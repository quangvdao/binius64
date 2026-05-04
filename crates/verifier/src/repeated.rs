// Copyright 2026 The Binius Developers

//! Descriptor for verifier-side batching of repeated identical constraint systems.

use binius_core::constraint_system::{
	AndConstraint, ConstraintSystem, MulConstraint, Operand, ShiftedValueIndex, ValueIndex,
	ValueVecLayout,
};
use binius_utils::SerializeBytes;
use bytes::BytesMut;
use digest::Digest;
use sha2::Sha256;

use crate::config::B128;

const REPEATED_DESCRIPTOR_VERSION: u128 = 1;
const REPEATED_DESCRIPTOR_HASH_TAG: &[u8] = b"binius64.repeated-constraint-system.v1";

/// Value-vector layout for repeated copies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RepeatedValueLayout {
	/// Instance-major tensor layout:
	///
	/// ```text
	/// flat_value = instance * base_value_count + local_value
	/// ```
	InstanceMajor = 0,
}

/// How per-instance public data is bound to the repeated proof.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RepeatedPublicBinding {
	/// Public inputs are still supplied through the flat public input section.
	FlatPublicInputs = 0,
}

/// A repeated identical constraint-system descriptor.
#[derive(Debug, Clone)]
pub struct RepeatedConstraintSystem {
	base: ConstraintSystem,
	log_instances: usize,
	value_layout: RepeatedValueLayout,
	public_binding: RepeatedPublicBinding,
	base_circuit_digest: [u8; 32],
}

impl RepeatedConstraintSystem {
	/// Constructs a descriptor using the v0 instance-major layout and flat public inputs.
	pub fn new(base: ConstraintSystem, log_instances: usize) -> Self {
		Self::with_modes(
			base,
			log_instances,
			RepeatedValueLayout::InstanceMajor,
			RepeatedPublicBinding::FlatPublicInputs,
		)
	}

	/// Constructs a descriptor with explicit layout and public binding modes.
	pub fn with_modes(
		base: ConstraintSystem,
		log_instances: usize,
		value_layout: RepeatedValueLayout,
		public_binding: RepeatedPublicBinding,
	) -> Self {
		assert!(base.value_vec_layout.committed_total_len.is_power_of_two());
		assert!(base.and_constraints.len().is_power_of_two());
		assert!(base.mul_constraints.len().is_power_of_two());
		let base_circuit_digest = digest_constraint_system(&base);
		Self {
			base,
			log_instances,
			value_layout,
			public_binding,
			base_circuit_digest,
		}
	}

	/// Returns the base constraint system.
	pub fn base(&self) -> &ConstraintSystem {
		&self.base
	}

	/// Returns the number of instance-axis variables.
	pub fn log_instances(&self) -> usize {
		self.log_instances
	}

	/// Returns the value layout mode.
	pub fn value_layout(&self) -> RepeatedValueLayout {
		self.value_layout
	}

	/// Returns the public binding mode.
	pub fn public_binding(&self) -> RepeatedPublicBinding {
		self.public_binding
	}

	/// Returns the SHA-256 digest of the serialized base constraint system.
	pub fn base_circuit_digest(&self) -> [u8; 32] {
		self.base_circuit_digest
	}

	/// Returns transcript elements binding the repeated descriptor.
	pub fn binding_scalars(&self) -> [B128; 8] {
		let [digest_lo, digest_hi] = digest_to_b128s(self.base_circuit_digest);
		[
			B128::new(REPEATED_DESCRIPTOR_VERSION),
			B128::new(self.log_instances as u128),
			B128::new(self.value_layout as u128),
			B128::new(self.public_binding as u128),
			B128::new(self.base.value_vec_layout.committed_total_len as u128),
			B128::new(
				(self.base.and_constraints.len() as u128)
					| ((self.base.mul_constraints.len() as u128) << 64),
			),
			digest_lo,
			digest_hi,
		]
	}

	/// Returns whether `flat` has the tensor dimensions implied by this repeated descriptor.
	pub fn matches_flat_shape(&self, flat: &ConstraintSystem) -> bool {
		let instances = 1usize << self.log_instances;
		flat.value_vec_layout.committed_total_len
			== self.base.value_vec_layout.committed_total_len * instances
			&& flat.and_constraints.len() == self.base.and_constraints.len() * instances
			&& flat.mul_constraints.len() == self.base.mul_constraints.len() * instances
	}

	/// Expands the repeated descriptor into the flat constraint system represented by v0.
	///
	/// The v0 instance-major layout repeats committed value indices by adding
	/// `instance * base_value_count` to every term. It keeps only the base public prefix public in
	/// the flat circuit, so values in later instances are hidden witness data. Base constants are
	/// therefore not supported by this layout: offsetting a constant for later instances would
	/// turn it into unconstrained private data.
	pub fn to_flat_constraint_system(&self) -> ConstraintSystem {
		assert_eq!(self.value_layout, RepeatedValueLayout::InstanceMajor);
		assert_eq!(self.public_binding, RepeatedPublicBinding::FlatPublicInputs);
		assert!(
			self.base.constants.is_empty(),
			"v0 repeated expansion requires circuits with no base constants"
		);

		let base_value_count = self.base.value_vec_layout.committed_total_len;
		ConstraintSystem::new(
			self.base.constants.clone(),
			self.repeated_value_vec_layout(),
			repeat_and_constraints(
				&self.base.and_constraints,
				self.log_instances,
				base_value_count,
			),
			repeat_mul_constraints(
				&self.base.mul_constraints,
				self.log_instances,
				base_value_count,
			),
		)
	}

	/// Returns whether `flat` is exactly the v0 flat expansion of this repeated descriptor.
	///
	/// This is intended for setup-time guardrails. It is linear in the flat constraint-system size,
	/// so hot verifier paths should bind the descriptor and then use the structured repeated
	/// verifier rather than recomputing this check every proof.
	pub fn matches_flat_constraint_system(&self, flat: &ConstraintSystem) -> bool {
		if self.value_layout != RepeatedValueLayout::InstanceMajor
			|| self.public_binding != RepeatedPublicBinding::FlatPublicInputs
			|| !self.base.constants.is_empty()
			|| !self.matches_flat_shape(flat)
			|| flat.constants != self.base.constants
			|| flat.value_vec_layout != self.repeated_value_vec_layout()
		{
			return false;
		}

		let base_value_count = self.base.value_vec_layout.committed_total_len;
		repeated_and_constraints_match(
			&self.base.and_constraints,
			&flat.and_constraints,
			self.log_instances,
			base_value_count,
		) && repeated_mul_constraints_match(
			&self.base.mul_constraints,
			&flat.mul_constraints,
			self.log_instances,
			base_value_count,
		)
	}

	fn repeated_value_vec_layout(&self) -> ValueVecLayout {
		let base_layout = &self.base.value_vec_layout;
		let repeated_total_len = base_layout.committed_total_len << self.log_instances;
		ValueVecLayout {
			n_const: base_layout.n_const,
			n_inout: base_layout.n_inout,
			n_witness: repeated_total_len - base_layout.offset_witness,
			n_internal: 0,
			offset_inout: base_layout.offset_inout,
			offset_witness: base_layout.offset_witness,
			committed_total_len: repeated_total_len,
			n_scratch: base_layout.n_scratch << self.log_instances,
		}
	}
}

fn repeat_and_constraints(
	base_constraints: &[AndConstraint],
	log_instances: usize,
	base_value_count: usize,
) -> Vec<AndConstraint> {
	let instances = 1 << log_instances;
	let mut constraints = Vec::with_capacity(instances * base_constraints.len());
	for instance in 0..instances {
		let value_offset = instance * base_value_count;
		for constraint in base_constraints {
			constraints.push(AndConstraint {
				a: offset_operand(&constraint.a, value_offset),
				b: offset_operand(&constraint.b, value_offset),
				c: offset_operand(&constraint.c, value_offset),
			});
		}
	}
	constraints
}

fn repeat_mul_constraints(
	base_constraints: &[MulConstraint],
	log_instances: usize,
	base_value_count: usize,
) -> Vec<MulConstraint> {
	let instances = 1 << log_instances;
	let mut constraints = Vec::with_capacity(instances * base_constraints.len());
	for instance in 0..instances {
		let value_offset = instance * base_value_count;
		for constraint in base_constraints {
			constraints.push(MulConstraint {
				a: offset_operand(&constraint.a, value_offset),
				b: offset_operand(&constraint.b, value_offset),
				lo: offset_operand(&constraint.lo, value_offset),
				hi: offset_operand(&constraint.hi, value_offset),
			});
		}
	}
	constraints
}

fn repeated_and_constraints_match(
	base_constraints: &[AndConstraint],
	flat_constraints: &[AndConstraint],
	log_instances: usize,
	base_value_count: usize,
) -> bool {
	let base_constraint_count = base_constraints.len();
	let instances = 1 << log_instances;
	if flat_constraints.len() != base_constraint_count * instances {
		return false;
	}
	for instance in 0..instances {
		let value_offset = instance * base_value_count;
		for (base_index, base_constraint) in base_constraints.iter().enumerate() {
			let flat_constraint = &flat_constraints[instance * base_constraint_count + base_index];
			if !offset_operand_matches(&base_constraint.a, &flat_constraint.a, value_offset)
				|| !offset_operand_matches(&base_constraint.b, &flat_constraint.b, value_offset)
				|| !offset_operand_matches(&base_constraint.c, &flat_constraint.c, value_offset)
			{
				return false;
			}
		}
	}
	true
}

fn repeated_mul_constraints_match(
	base_constraints: &[MulConstraint],
	flat_constraints: &[MulConstraint],
	log_instances: usize,
	base_value_count: usize,
) -> bool {
	let base_constraint_count = base_constraints.len();
	let instances = 1 << log_instances;
	if flat_constraints.len() != base_constraint_count * instances {
		return false;
	}
	for instance in 0..instances {
		let value_offset = instance * base_value_count;
		for (base_index, base_constraint) in base_constraints.iter().enumerate() {
			let flat_constraint = &flat_constraints[instance * base_constraint_count + base_index];
			if !offset_operand_matches(&base_constraint.a, &flat_constraint.a, value_offset)
				|| !offset_operand_matches(&base_constraint.b, &flat_constraint.b, value_offset)
				|| !offset_operand_matches(&base_constraint.lo, &flat_constraint.lo, value_offset)
				|| !offset_operand_matches(&base_constraint.hi, &flat_constraint.hi, value_offset)
			{
				return false;
			}
		}
	}
	true
}

fn offset_operand(operand: &Operand, offset: usize) -> Operand {
	operand
		.iter()
		.map(|term| offset_term(term, offset))
		.collect()
}

fn offset_operand_matches(base: &Operand, flat: &Operand, offset: usize) -> bool {
	base.len() == flat.len()
		&& std::iter::zip(base, flat).all(|(base_term, flat_term)| {
			let expected = offset_term(base_term, offset);
			expected.value_index == flat_term.value_index
				&& expected.shift_variant == flat_term.shift_variant
				&& expected.amount == flat_term.amount
		})
}

fn offset_term(term: &ShiftedValueIndex, offset: usize) -> ShiftedValueIndex {
	ShiftedValueIndex {
		value_index: ValueIndex(term.value_index.0 + offset as u32),
		shift_variant: term.shift_variant,
		amount: term.amount,
	}
}

fn digest_constraint_system(constraint_system: &ConstraintSystem) -> [u8; 32] {
	let mut bytes = BytesMut::new();
	constraint_system
		.serialize(&mut bytes)
		.expect("serializing constraint system into bytes should not fail");

	let mut hasher = Sha256::new();
	hasher.update(REPEATED_DESCRIPTOR_HASH_TAG);
	hasher.update(&bytes);
	hasher.finalize().into()
}

fn digest_to_b128s(digest: [u8; 32]) -> [B128; 2] {
	let mut lo = [0u8; 16];
	let mut hi = [0u8; 16];
	lo.copy_from_slice(&digest[..16]);
	hi.copy_from_slice(&digest[16..]);
	[
		B128::new(u128::from_le_bytes(lo)),
		B128::new(u128::from_le_bytes(hi)),
	]
}

#[cfg(test)]
mod tests {
	use binius_core::{
		constraint_system::{AndConstraint, ConstraintSystem, MulConstraint, ValueVecLayout},
		word::Word,
	};

	use super::*;

	fn value_vec_layout(value_count: usize) -> ValueVecLayout {
		ValueVecLayout {
			n_const: 0,
			n_inout: 2,
			n_witness: value_count - 2,
			n_internal: 0,
			offset_inout: 0,
			offset_witness: 2,
			committed_total_len: value_count,
			n_scratch: 0,
		}
	}

	fn base_constraint_system() -> ConstraintSystem {
		ConstraintSystem::new(
			Vec::new(),
			value_vec_layout(4),
			vec![AndConstraint::default()],
			vec![MulConstraint::default()],
		)
	}

	#[test]
	fn flat_expansion_matches_repeated_descriptor() {
		let repeated = RepeatedConstraintSystem::new(base_constraint_system(), 2);
		let flat = repeated.to_flat_constraint_system();

		assert!(repeated.matches_flat_shape(&flat));
		assert!(repeated.matches_flat_constraint_system(&flat));
		assert_eq!(flat.value_vec_layout.committed_total_len, 16);
		assert_eq!(flat.value_vec_layout.offset_witness, 2);
		assert_eq!(flat.and_constraints.len(), 4);
		assert_eq!(flat.mul_constraints.len(), 4);
	}

	#[test]
	fn flat_expansion_rejects_same_shape_different_structure() {
		let repeated = RepeatedConstraintSystem::new(base_constraint_system(), 2);
		let mut flat = repeated.to_flat_constraint_system();
		flat.and_constraints[1]
			.a
			.push(ShiftedValueIndex::plain(ValueIndex(0)));

		assert!(repeated.matches_flat_shape(&flat));
		assert!(!repeated.matches_flat_constraint_system(&flat));
	}

	#[test]
	#[should_panic(expected = "no base constants")]
	fn flat_expansion_rejects_base_constants() {
		let repeated = RepeatedConstraintSystem::new(
			ConstraintSystem::new(
				vec![Word::from_u64(1)],
				ValueVecLayout {
					n_const: 1,
					n_inout: 1,
					n_witness: 2,
					n_internal: 0,
					offset_inout: 1,
					offset_witness: 2,
					committed_total_len: 4,
					n_scratch: 0,
				},
				vec![AndConstraint::default()],
				vec![MulConstraint::default()],
			),
			1,
		);
		let _ = repeated.to_flat_constraint_system();
	}
}
