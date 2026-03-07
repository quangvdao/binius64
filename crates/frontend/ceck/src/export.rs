// Copyright 2026 The Binius Developers
//! Exports a `ConstraintSystem` to `.ceck` S-expression format.
//!
//! This enables formal equivalence checking between compiled circuits
//! and hand-written reference specifications.

use binius_core::constraint_system::{ConstraintSystem, Operand, ShiftVariant, ShiftedValueIndex};

/// Export a `ConstraintSystem` to a `.ceck` S-expression string.
///
/// Wire names follow the convention:
/// - Constants: `0x{value:016x}`
/// - Witnesses: `$w{i}` where `i` is the witness index (0-based)
pub fn export_constraint_set(cs: &ConstraintSystem) -> String {
	let mut out = String::from("(constraint_set\n");

	for c in &cs.and_constraints {
		out.push_str("  (and ");
		out.push_str(&format_operand(cs, &c.a));
		out.push(' ');
		out.push_str(&format_operand(cs, &c.b));
		out.push(' ');
		out.push_str(&format_operand(cs, &c.c));
		out.push_str(")\n");
	}

	for c in &cs.mul_constraints {
		out.push_str("  (mul ");
		out.push_str(&format_operand(cs, &c.a));
		out.push(' ');
		out.push_str(&format_operand(cs, &c.b));
		out.push(' ');
		out.push_str(&format_operand(cs, &c.hi));
		out.push(' ');
		out.push_str(&format_operand(cs, &c.lo));
		out.push_str(")\n");
	}

	out.push(')');
	out
}

/// Export two constraint systems wrapped in an `(assert_eqv ...)` block.
pub fn export_assert_eqv(lhs: &ConstraintSystem, rhs: &ConstraintSystem) -> String {
	format!(
		"(assert_eqv\n  {}\n  {}\n)",
		export_constraint_set(lhs),
		export_constraint_set(rhs),
	)
}

fn format_operand(cs: &ConstraintSystem, operand: &Operand) -> String {
	if operand.len() == 1 {
		format_shifted_value(cs, &operand[0])
	} else {
		let terms: Vec<String> = operand.iter().map(|sv| format_shifted_value(cs, sv)).collect();
		format!("(xor {})", terms.join(" "))
	}
}

fn format_shifted_value(cs: &ConstraintSystem, sv: &ShiftedValueIndex) -> String {
	let base = format_value_index(cs, sv.value_index.0 as usize);

	if sv.amount == 0 {
		return base;
	}

	let op = match sv.shift_variant {
		ShiftVariant::Sll => "sll",
		ShiftVariant::Slr => "slr",
		ShiftVariant::Sar => "sar",
		ShiftVariant::Rotr => "ror",
		ShiftVariant::Sll32 => "sll32",
		ShiftVariant::Srl32 => "slr32",
		ShiftVariant::Sra32 => "sar32",
		ShiftVariant::Rotr32 => "ror32",
	};

	format!("({op} {base} {})", sv.amount)
}

fn format_value_index(cs: &ConstraintSystem, idx: usize) -> String {
	let n_const = cs.constants.len();
	let n_inout = cs.value_vec_layout.n_inout;
	let witness_start = cs.value_vec_layout.offset_witness;

	if idx < n_const {
		format!("0x{:016x}", cs.constants[idx].0)
	} else if idx < n_const + n_inout {
		format!("$io{}", idx - n_const)
	} else if idx >= witness_start {
		format!("$w{}", idx - witness_start)
	} else {
		format!("$internal{idx}")
	}
}

#[cfg(test)]
mod tests {
	use binius_core::{
		constraint_system::{AndConstraint, ValueIndex, ValueVecLayout},
		word::Word,
	};

	use super::*;

	#[test]
	fn test_export_simple_and() {
		let layout = ValueVecLayout {
			n_const: 1,
			n_inout: 0,
			n_witness: 3,
			n_internal: 0,
			offset_inout: 1,
			offset_witness: 2,
			committed_total_len: 8,
			n_scratch: 0,
		};

		let constraint = AndConstraint::plain_abc(
			vec![ValueIndex(2)],
			vec![ValueIndex(3)],
			vec![ValueIndex(4)],
		);

		let cs = ConstraintSystem::new(vec![Word(0)], layout, vec![constraint], vec![]);

		let exported = export_constraint_set(&cs);
		assert!(exported.contains("(and $w0 $w1 $w2)"));
	}

	#[test]
	fn test_export_with_shifts() {
		let layout = ValueVecLayout {
			n_const: 1,
			n_inout: 0,
			n_witness: 2,
			n_internal: 0,
			offset_inout: 1,
			offset_witness: 2,
			committed_total_len: 4,
			n_scratch: 0,
		};

		let constraint = AndConstraint {
			a: vec![ShiftedValueIndex::rotr32(ValueIndex(2), 16)],
			b: vec![ShiftedValueIndex::plain(ValueIndex(0))],
			c: vec![ShiftedValueIndex::plain(ValueIndex(3))],
		};

		let cs = ConstraintSystem::new(vec![Word(u64::MAX)], layout, vec![constraint], vec![]);

		let exported = export_constraint_set(&cs);
		assert!(exported.contains("(ror32 $w0 16)"), "exported: {exported}");
		assert!(
			exported.contains("0xffffffffffffffff"),
			"exported: {exported}"
		);
	}
}
