import pytest
import json
import tempfile
from pathlib import Path
import main 

def test_dump_genome(capsys):
    """test that dump_genome prints JSON matching an expected file."""

    gen_ref = main.build_genome_reference("test/main_test/smpref_T1_5uniq_genome.fa", 31)
    # capture the printed output of dump_genome func
    main.dump_genome(gen_ref)
    captured = capsys.readouterr()  # get output
    # Convert the captured output to JSON format
    output_json = json.loads(captured.out.strip())
    with open("test/main_test/smpref_T1_5uniq_genome_dumpref.json", "r") as f:
        expected_json = json.load(f)
    assert output_json == expected_json

    gen_ref = main.build_genome_reference("test/main_test/smpref_T2_3uniq_2ambig_genome.fa", 31)
    main.dump_genome(gen_ref)
    captured = capsys.readouterr()
    output_json = json.loads(captured.out.strip())
    with open("test/main_test/smpref_T2_3uniq_2ambig_genome_dumpref.json", "r") as f:
        expected_json = json.load(f)
    assert output_json == expected_json

    gen_ref = main.build_genome_reference("test/main_test/smpref_T3_5ambig_genome.fa", 31)
    main.dump_genome(gen_ref)
    captured = capsys.readouterr()
    output_json = json.loads(captured.out.strip())
    with open("test/main_test/smpref_T3_5ambig_genome_dumpref.json", "r") as f:
        expected_json = json.load(f)
    assert output_json == expected_json

    gen_ref = main.build_genome_reference("test/main_test/smpref_T4_7uniq_3ambig_genome.fa", 31)
    main.dump_genome(gen_ref)
    captured = capsys.readouterr()
    output_json = json.loads(captured.out.strip())
    with open("test/main_test/smpref_T4_7uniq_3ambig_genome_dumpref.json", "r") as f:
        expected_json = json.load(f)
    assert output_json == expected_json

    gen_ref = main.build_genome_reference("test/main_test/smpref_T5_5uniq_1k_genome.fa", 31)
    main.dump_genome(gen_ref)
    captured = capsys.readouterr()
    output_json = json.loads(captured.out.strip())
    with open("test/main_test/smpref_T5_5uniq_1k_genome_dumpref.json", "r") as f:
        expected_json = json.load(f)
    assert output_json == expected_json


def test_dump_align(capsys):
    """Test that dump_align prints JSON matching an expected file."""
    gen_ref = main.build_genome_reference("example_data/example_ref_mid_clean.fa", 150)
    align = main.align("example_data/mid_example_se_k150_n10k.fq", gen_ref, 1, 5, None)
    # Capture printed output
    main.dump_align(align)
    captured = capsys.readouterr()
    # Convert captured output to JSON format
    output_json = json.loads(captured.out)
    # Load expected JSON
    with open("test/main_test/mid_example_se_k150_n10k_dump.json", "r") as f:
        expected_json = json.load(f)
    assert output_json == expected_json
