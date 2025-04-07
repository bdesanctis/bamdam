#!/usr/bin/env python3
"""
Simple tests for bamdam commands.
"""

import pytest
import argparse
from pathlib import Path


# Import the command functions directly from bamdam
from bamdam.bamdam import (
    shrink,
    compute,
    extract,
    plotdamage,
    plotbaminfo,
    combine,
    krona,
)


def test_shrink(tmp_path):
    """Test the shrink command."""
    # Download test data if needed
    sample_lca = "tests/data/small.lca"
    sample_bam = "tests/data/small.bam"

    # Set up output paths in temp directory
    out_lca = tmp_path / "small.lca"
    out_bam = tmp_path / "small.bam"

    # Create args and call function
    args = argparse.Namespace()
    args.in_lca = str(sample_lca)
    args.in_bam = str(sample_bam)
    args.out_lca = str(out_lca)
    args.out_bam = str(out_bam)
    args.stranded = "ds"
    args.upto = "family"
    args.mincount = 1
    args.minsim = 0.05
    args.exclude_keywords = []
    args.exclude_keyword_file = None
    args.annotate_pmd = False

    shrink(args)

    # Check if output files exist
    assert out_lca.exists()
    assert out_bam.exists()
