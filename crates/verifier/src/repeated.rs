// Copyright 2026 The Binius Developers

//! Descriptor for verifier-side batching of repeated identical constraint systems.

use binius_core::constraint_system::ConstraintSystem;
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
