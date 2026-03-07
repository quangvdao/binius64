// Copyright 2026 The Binius Developers

use anyhow::Result;
use binius_circuits::blake3::Blake3;
use binius_frontend::{CircuitBuilder, WitnessFiller};
use clap::Args;

use super::utils;
use crate::ExampleCircuit;

/// BLAKE3 circuit example demonstrating the BLAKE3 hash function implementation.
///
/// Supports batching multiple independent hash invocations into a single circuit
/// via the `count` parameter for amortized proving cost.
pub struct Blake3Example {
	gadgets: Vec<Blake3>,
}

/// Circuit parameters that affect structure (compile-time configuration)
#[derive(Debug, Clone, Args)]
pub struct Params {
	/// Maximum message length in bytes per invocation. Messages over 1024 bytes
	/// use tree hashing with multiple chunks.
	#[arg(long)]
	pub max_bytes: Option<usize>,

	/// Number of independent BLAKE3 hash invocations to batch.
	#[arg(long, default_value = "1")]
	pub count: usize,
}

/// Instance data for witness population (runtime values)
#[derive(Debug, Clone, Args)]
#[group(multiple = false)]
pub struct Instance {
	/// Length of the randomly generated message, in bytes (defaults to 1024).
	#[arg(long)]
	pub message_len: Option<usize>,

	/// UTF-8 string to hash (if not provided, random bytes are generated)
	#[arg(long)]
	pub message_string: Option<String>,
}

impl ExampleCircuit for Blake3Example {
	type Params = Params;
	type Instance = Instance;

	fn build(params: Params, builder: &mut CircuitBuilder) -> Result<Self> {
		let max_bytes = utils::determine_hash_max_bytes_from_args(params.max_bytes)?;

		let gadgets = (0..params.count)
			.map(|_| Blake3::new_witness(builder, max_bytes))
			.collect();

		Ok(Self { gadgets })
	}

	fn populate_witness(&self, instance: Instance, w: &mut WitnessFiller) -> Result<()> {
		let msg_len = self.gadgets[0].length;

		for gadget in &self.gadgets {
			let raw_message =
				utils::generate_message_bytes(instance.message_string.clone(), instance.message_len);
			let padded_message = utils::zero_pad_message(raw_message, msg_len)?;

			let digest: [u8; 32] = blake3::hash(&padded_message).into();

			gadget.populate_message(w, &padded_message);
			gadget.populate_digest(w, &digest);
		}

		Ok(())
	}

	fn param_summary(params: &Self::Params) -> Option<String> {
		let bytes = params
			.max_bytes
			.unwrap_or(utils::DEFAULT_HASH_MESSAGE_BYTES);
		if params.count > 1 {
			Some(format!("{}x{}b", params.count, bytes))
		} else {
			Some(format!("{bytes}b"))
		}
	}
}
