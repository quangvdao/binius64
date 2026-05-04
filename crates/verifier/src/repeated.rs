// Copyright 2026 The Binius Developers

//! Descriptor for verifier-side batching of repeated identical constraint systems.

use binius_core::constraint_system::{
	AndConstraint, ConstraintSystem, MulConstraint, Operand, ShiftedValueIndex, ValueIndex,
	ValueVec, ValueVecLayout,
};
use binius_core::error::ConstraintSystemError;
use binius_core::word::Word;
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
	/// flat_value = local_value                                  if local_value is a constant
	/// flat_value = instance * base_value_count + local_value    otherwise
	/// ```
	///
	/// Base constants are shared across every repeated instance. All non-constant committed
	/// values, including per-instance public inputs in the base circuit, are repeated
	/// instance-major. With [`RepeatedPublicBinding::FlatPublicInputs`], only the first
	/// instance's public prefix remains public in the flat circuit; later instance values are
	/// hidden witness words.
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

/// Error returned when flattening per-instance value vectors for a repeated descriptor.
#[derive(Debug, thiserror::Error)]
pub enum RepeatedValueVecError {
	#[error("incorrect repeated instance count: expected {expected}, got {actual}")]
	IncorrectInstanceCount { expected: usize, actual: usize },
	#[error(
		"instance {instance} committed value length mismatch: expected {expected}, got {actual}"
	)]
	InstanceValueVecLen {
		instance: usize,
		expected: usize,
		actual: usize,
	},
	#[error("instance {instance} public value length mismatch: expected {expected}, got {actual}")]
	InstancePublicLen {
		instance: usize,
		expected: usize,
		actual: usize,
	},
	#[error("instance {instance} has a different shared constant at index {index}")]
	SharedConstantMismatch { instance: usize, index: usize },
	#[error("value-vector construction error: {0}")]
	ValueVec(#[from] ConstraintSystemError),
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
	/// The v0 instance-major layout shares base constants across every repeated instance and
	/// repeats all non-constant committed value indices by adding
	/// `instance * base_value_count` to every term. It keeps only the base public prefix public in
	/// the flat circuit, so values in later instances are hidden witness data.
	pub fn to_flat_constraint_system(&self) -> ConstraintSystem {
		assert_eq!(self.value_layout, RepeatedValueLayout::InstanceMajor);
		assert_eq!(self.public_binding, RepeatedPublicBinding::FlatPublicInputs);

		let base_value_count = self.base.value_vec_layout.committed_total_len;
		let base_const_count = self.base.value_vec_layout.n_const;
		ConstraintSystem::new(
			self.base.constants.clone(),
			self.repeated_value_vec_layout(),
			repeat_and_constraints(
				&self.base.and_constraints,
				self.log_instances,
				base_value_count,
				base_const_count,
			),
			repeat_mul_constraints(
				&self.base.mul_constraints,
				self.log_instances,
				base_value_count,
				base_const_count,
			),
		)
	}

	/// Flattens one base-layout value vector per instance into the v0 flat repeated layout.
	///
	/// Constants are shared and must match the first instance. All non-constant committed values
	/// are copied into the instance-major slot used by [`Self::to_flat_constraint_system`]. The
	/// resulting value vector exposes only the flat circuit public prefix; every later-instance
	/// value lives in the non-public portion.
	pub fn to_flat_value_vec(
		&self,
		instances: &[ValueVec],
	) -> Result<ValueVec, RepeatedValueVecError> {
		assert_eq!(self.value_layout, RepeatedValueLayout::InstanceMajor);
		assert_eq!(self.public_binding, RepeatedPublicBinding::FlatPublicInputs);

		let expected_instances = 1usize << self.log_instances;
		if instances.len() != expected_instances {
			return Err(RepeatedValueVecError::IncorrectInstanceCount {
				expected: expected_instances,
				actual: instances.len(),
			});
		}

		let base_layout = &self.base.value_vec_layout;
		let base_value_count = base_layout.committed_total_len;
		for (instance, value_vec) in instances.iter().enumerate() {
			if value_vec.size() != base_value_count {
				return Err(RepeatedValueVecError::InstanceValueVecLen {
					instance,
					expected: base_value_count,
					actual: value_vec.size(),
				});
			}
			if value_vec.public().len() != base_layout.offset_witness {
				return Err(RepeatedValueVecError::InstancePublicLen {
					instance,
					expected: base_layout.offset_witness,
					actual: value_vec.public().len(),
				});
			}
			for constant_index in 0..base_layout.n_const {
				if value_vec.get(constant_index) != instances[0].get(constant_index) {
					return Err(RepeatedValueVecError::SharedConstantMismatch {
						instance,
						index: constant_index,
					});
				}
			}
		}

		let flat_layout = self.repeated_value_vec_layout();
		let mut flat_values = vec![Word::ZERO; flat_layout.committed_total_len];
		for (instance_index, value_vec) in instances.iter().enumerate() {
			for local_value_index in 0..base_value_count {
				let flat_value_index = if local_value_index < base_layout.n_const {
					local_value_index
				} else {
					instance_index * base_value_count + local_value_index
				};
				flat_values[flat_value_index] = value_vec.get(local_value_index);
			}
		}

		let private = flat_values.split_off(flat_layout.offset_witness);
		ValueVec::new_from_data(flat_layout, flat_values, private).map_err(Into::into)
	}

	/// Returns whether `flat` is exactly the v0 flat expansion of this repeated descriptor.
	///
	/// This is intended for setup-time guardrails. It is linear in the flat constraint-system size,
	/// so hot verifier paths should bind the descriptor and then use the structured repeated
	/// verifier rather than recomputing this check every proof.
	pub fn matches_flat_constraint_system(&self, flat: &ConstraintSystem) -> bool {
		if self.value_layout != RepeatedValueLayout::InstanceMajor
			|| self.public_binding != RepeatedPublicBinding::FlatPublicInputs
			|| !self.matches_flat_shape(flat)
			|| flat.constants != self.base.constants
			|| flat.value_vec_layout != self.repeated_value_vec_layout()
		{
			return false;
		}

		let base_value_count = self.base.value_vec_layout.committed_total_len;
		let base_const_count = self.base.value_vec_layout.n_const;
		repeated_and_constraints_match(
			&self.base.and_constraints,
			&flat.and_constraints,
			self.log_instances,
			base_value_count,
			base_const_count,
		) && repeated_mul_constraints_match(
			&self.base.mul_constraints,
			&flat.mul_constraints,
			self.log_instances,
			base_value_count,
			base_const_count,
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
	base_const_count: usize,
) -> Vec<AndConstraint> {
	let instances = 1 << log_instances;
	let mut constraints = Vec::with_capacity(instances * base_constraints.len());
	for instance in 0..instances {
		let value_offset = instance * base_value_count;
		for constraint in base_constraints {
			constraints.push(AndConstraint {
				a: offset_operand(&constraint.a, value_offset, base_const_count),
				b: offset_operand(&constraint.b, value_offset, base_const_count),
				c: offset_operand(&constraint.c, value_offset, base_const_count),
			});
		}
	}
	constraints
}

fn repeat_mul_constraints(
	base_constraints: &[MulConstraint],
	log_instances: usize,
	base_value_count: usize,
	base_const_count: usize,
) -> Vec<MulConstraint> {
	let instances = 1 << log_instances;
	let mut constraints = Vec::with_capacity(instances * base_constraints.len());
	for instance in 0..instances {
		let value_offset = instance * base_value_count;
		for constraint in base_constraints {
			constraints.push(MulConstraint {
				a: offset_operand(&constraint.a, value_offset, base_const_count),
				b: offset_operand(&constraint.b, value_offset, base_const_count),
				lo: offset_operand(&constraint.lo, value_offset, base_const_count),
				hi: offset_operand(&constraint.hi, value_offset, base_const_count),
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
	base_const_count: usize,
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
			if !offset_operand_matches(
				&base_constraint.a,
				&flat_constraint.a,
				value_offset,
				base_const_count,
			) || !offset_operand_matches(
				&base_constraint.b,
				&flat_constraint.b,
				value_offset,
				base_const_count,
			) || !offset_operand_matches(
				&base_constraint.c,
				&flat_constraint.c,
				value_offset,
				base_const_count,
			) {
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
	base_const_count: usize,
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
			if !offset_operand_matches(
				&base_constraint.a,
				&flat_constraint.a,
				value_offset,
				base_const_count,
			) || !offset_operand_matches(
				&base_constraint.b,
				&flat_constraint.b,
				value_offset,
				base_const_count,
			) || !offset_operand_matches(
				&base_constraint.lo,
				&flat_constraint.lo,
				value_offset,
				base_const_count,
			) || !offset_operand_matches(
				&base_constraint.hi,
				&flat_constraint.hi,
				value_offset,
				base_const_count,
			) {
				return false;
			}
		}
	}
	true
}

fn offset_operand(operand: &Operand, offset: usize, base_const_count: usize) -> Operand {
	operand
		.iter()
		.map(|term| offset_term(term, offset, base_const_count))
		.collect()
}

fn offset_operand_matches(
	base: &Operand,
	flat: &Operand,
	offset: usize,
	base_const_count: usize,
) -> bool {
	base.len() == flat.len()
		&& std::iter::zip(base, flat).all(|(base_term, flat_term)| {
			let expected = offset_term(base_term, offset, base_const_count);
			expected.value_index == flat_term.value_index
				&& expected.shift_variant == flat_term.shift_variant
				&& expected.amount == flat_term.amount
		})
}

fn offset_term(
	term: &ShiftedValueIndex,
	offset: usize,
	base_const_count: usize,
) -> ShiftedValueIndex {
	let value_index = if term.value_index.0 < base_const_count as u32 {
		term.value_index
	} else {
		ValueIndex(term.value_index.0 + offset as u32)
	};
	ShiftedValueIndex {
		value_index,
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
		constraint_system::{
			AndConstraint, ConstraintSystem, MulConstraint, ValueVec, ValueVecLayout,
		},
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

	fn base_constraint_system_with_constants() -> ConstraintSystem {
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
			vec![AndConstraint {
				a: vec![ShiftedValueIndex::plain(ValueIndex(0))],
				b: vec![ShiftedValueIndex::plain(ValueIndex(2))],
				c: vec![ShiftedValueIndex::plain(ValueIndex(3))],
			}],
			vec![MulConstraint::default()],
		)
	}

	fn base_value_vec_with_constant(inout: u64, witness_0: u64, witness_1: u64) -> ValueVec {
		ValueVec::new_from_data(
			base_constraint_system_with_constants().value_vec_layout,
			vec![Word::from_u64(1), Word::from_u64(inout)],
			vec![Word::from_u64(witness_0), Word::from_u64(witness_1)],
		)
		.expect("base value vec has matching layout")
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
	fn flat_expansion_shares_base_constants() {
		let repeated = RepeatedConstraintSystem::new(base_constraint_system_with_constants(), 2);
		let flat = repeated.to_flat_constraint_system();

		assert!(repeated.matches_flat_shape(&flat));
		assert!(repeated.matches_flat_constraint_system(&flat));
		assert_eq!(flat.constants, vec![Word::from_u64(1)]);
		assert_eq!(flat.value_vec_layout.n_const, 1);
		assert_eq!(flat.and_constraints[0].a[0].value_index, ValueIndex(0));
		assert_eq!(flat.and_constraints[1].a[0].value_index, ValueIndex(0));
		assert_eq!(flat.and_constraints[1].b[0].value_index, ValueIndex(6));
		assert_eq!(flat.and_constraints[1].c[0].value_index, ValueIndex(7));
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
	fn flat_expansion_rejects_duplicated_constants() {
		let repeated = RepeatedConstraintSystem::new(base_constraint_system_with_constants(), 1);
		let mut flat = repeated.to_flat_constraint_system();
		flat.and_constraints[1].a[0].value_index = ValueIndex(4);

		assert!(repeated.matches_flat_shape(&flat));
		assert!(!repeated.matches_flat_constraint_system(&flat));
	}

	#[test]
	fn flat_value_vec_shares_constants_and_hides_later_instances() {
		let repeated = RepeatedConstraintSystem::new(base_constraint_system_with_constants(), 2);
		let instances = vec![
			base_value_vec_with_constant(10, 20, 30),
			base_value_vec_with_constant(11, 21, 31),
			base_value_vec_with_constant(12, 22, 32),
			base_value_vec_with_constant(13, 23, 33),
		];

		let flat = repeated
			.to_flat_value_vec(&instances)
			.expect("instances flatten");

		assert_eq!(flat.public(), &[Word::from_u64(1), Word::from_u64(10)]);
		assert_eq!(flat.get(0), Word::from_u64(1));
		assert_eq!(flat.get(1), Word::from_u64(10));
		assert_eq!(flat.get(2), Word::from_u64(20));
		assert_eq!(flat.get(3), Word::from_u64(30));
		assert_eq!(flat.get(4), Word::ZERO);
		assert_eq!(flat.get(5), Word::from_u64(11));
		assert_eq!(flat.get(6), Word::from_u64(21));
		assert_eq!(flat.get(7), Word::from_u64(31));
		assert_eq!(flat.get(12), Word::ZERO);
		assert_eq!(flat.get(13), Word::from_u64(13));
		assert_eq!(flat.get(14), Word::from_u64(23));
		assert_eq!(flat.get(15), Word::from_u64(33));
	}

	#[test]
	fn flat_value_vec_rejects_mismatched_shared_constants() {
		let repeated = RepeatedConstraintSystem::new(base_constraint_system_with_constants(), 1);
		let mut mismatched = base_value_vec_with_constant(11, 21, 31);
		mismatched.set(0, Word::from_u64(2));

		let err = repeated
			.to_flat_value_vec(&[base_value_vec_with_constant(10, 20, 30), mismatched])
			.expect_err("mismatched constant is rejected");
		assert!(matches!(
			err,
			RepeatedValueVecError::SharedConstantMismatch {
				instance: 1,
				index: 0
			}
		));
	}
}
