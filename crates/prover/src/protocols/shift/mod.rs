// Copyright 2025 Irreducible Inc.

use binius_verifier::protocols::shift::{BITAND_ARITY, INTMUL_ARITY, SHIFT_VARIANT_COUNT};

mod error;
mod key_collection;
mod monster;
mod phase_1;
mod phase_2;
mod prove;

pub use error::Error;
pub(crate) use key_collection::ShiftKeySource;
pub use key_collection::{KeyCollection, build_key_collection};
pub(crate) use prove::prove_with_key_source;
pub use prove::{OperatorData, PreparedOperatorData, prove};
