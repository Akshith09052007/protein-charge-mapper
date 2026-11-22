import unittest

from rna_predictor import can_pair, predict_structure, validate_sequence


class TestRNASequence(unittest.TestCase):
    def test_whitespace_and_dna_t_are_normalized(self):
        self.assertEqual(validate_sequence(" atg c "), "AUGC")

    def test_invalid_and_empty_sequences(self):
        for sequence in ("", "   ", "ACGX"):
            with self.subTest(sequence=sequence):
                with self.assertRaises(ValueError):
                    validate_sequence(sequence)

    def test_pair_types(self):
        self.assertTrue(can_pair("A", "U"))
        self.assertTrue(can_pair("G", "C"))
        self.assertTrue(can_pair("G", "U"))
        self.assertFalse(can_pair("G", "U", allow_wobble=False))
        self.assertFalse(can_pair("A", "C"))


class TestPrediction(unittest.TestCase):
    def test_output_lengths_and_dot_bracket(self):
        result = predict_structure("GGGAAACCC")
        self.assertEqual(result["length"], len(result["sequence"]))
        self.assertEqual(len(result["structure"]), result["length"])
        self.assertEqual(result["structure"].count("("), result["base_pair_count"])
        self.assertEqual(result["structure"].count(")"), result["base_pair_count"])

    def test_minimum_loop_constraint(self):
        result = predict_structure("AUAU", min_loop=3, allow_wobble=False)
        self.assertEqual(result["base_pair_count"], 0)

    def test_wobble_setting(self):
        enabled = predict_structure("GGGAAACCC", min_loop=0, allow_wobble=True)
        disabled = predict_structure("GGGAAACCC", min_loop=0, allow_wobble=False)
        self.assertGreaterEqual(enabled["base_pair_count"], disabled["base_pair_count"])

    def test_single_nucleotide(self):
        result = predict_structure("A")
        self.assertEqual(result["structure"], ".")
        self.assertEqual(result["base_pair_count"], 0)

    def test_invalid_options(self):
        with self.assertRaises(ValueError):
            predict_structure("ACGU", min_loop=-1)
        with self.assertRaises(ValueError):
            predict_structure("ACGU", min_loop=1.5)
        with self.assertRaises(ValueError):
            predict_structure("ACGU", allow_wobble="false")


if __name__ == "__main__":
    unittest.main()
