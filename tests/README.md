# Test Architecture

Each test genome assembly should go into its own directory under `tests/data` with the desired
name as the genome. This way, you can define multiple variations of the genome assembly, such quality variations
or assembly graphs for testing of multiple scenarios.

Each genome should have defined databases for testing against, expected locus match and expected serotype.
We may need to also define expected confidence score in the future in case we want to test scenarios where
we expect an "Untypeable" result.