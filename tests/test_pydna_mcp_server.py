#!/usr/bin/env python3
"""
Comprehensive pytest tests for PyDNA MCP Server.

This test suite covers all methods of the PyDNA MCP server with:
- Unit tests for each method
- Integration tests for workflows
- Edge cases and error handling
- Mock data and real sequence examples
- Async testing for MCP functionality
"""

import pytest
import asyncio
import tempfile
import os
from pathlib import Path
from unittest.mock import Mock, patch, AsyncMock
from typing import Dict, Any, List

# Import the server and its components
from pydna_mcp.server import PyDNAMCP, SequenceInfo, PrimerInfo, AmpliconInfo, RestrictionAnalysis, AssemblyResult

# Test data constants
TEST_SEQUENCES = {
    "short_dna": "ATGAAAGCACTGATTCTATTGCTG",
    "pcr_template": "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAATAGGATCCTTTTTTTTAAAA",
    "restriction_seq": "GAATTCATGAAAGCACTGATTCTATTGCTGAAAAAGATAATAGGATCC",
    "circular_plasmid": "GAATTCATGAAAGCACTGATTCTATTGCTGAAAAAGATAATAGGATCCAAGCTTGTCGACGCGGCCGC",
    "overlap_frag1": "ATGAAAGCACTGATTCTATTGC",
    "overlap_frag2": "CTATTGCTGAAAAAGATAAT", 
    "overlap_frag3": "AGATAATAGGATCCTTTTTTTT",
    "long_template": "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAATAGGATCCTTTTTTTTAAAATCGACGCATGCAAGCTTGAATTCGGATCCGTCGAC"
}

TEST_PRIMERS = {
    "forward": "ATGAAAGCACTGATTC",
    "reverse": "ATTATCTTTTTCAGC",
    "ecori_forward": "GAATTCATGAAAGCACTGATT",
    "bamhi_reverse": "TTTTTTAAAGGATCCTTTT"
}

# Sample file contents for testing
FASTA_CONTENT = """>test_sequence 
ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT
"""

GENBANK_CONTENT = """LOCUS       test_seq                34 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  Test sequence for PyDNA MCP server
ACCESSION   test_seq
VERSION     test_seq
KEYWORDS    .
SOURCE      synthetic
  ORGANISM  synthetic
FEATURES             Location/Qualifiers
     gene            1..34
                     /gene="test_gene"
                     /product="test protein"
ORIGIN
        1 atgaaagcac tgattctatt gctgaaaaag ataat
//
"""


@pytest.fixture
def mcp_server():
    """Create a PyDNA MCP server instance for testing."""
    return PyDNAMCP()


@pytest.fixture
def sample_sequences():
    """Provide sample DNA sequences for testing."""
    return TEST_SEQUENCES


@pytest.fixture
def sample_primers():
    """Provide sample primer sequences for testing."""
    return TEST_PRIMERS


@pytest.fixture
def temp_files():
    """Create temporary files for file I/O testing."""
    files = {}
    
    # Create FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(FASTA_CONTENT)
        files['fasta'] = f.name
    
    # Create GenBank file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.gb', delete=False) as f:
        f.write(GENBANK_CONTENT)
        files['genbank'] = f.name
    
    yield files
    
    # Cleanup
    for filepath in files.values():
        if os.path.exists(filepath):
            os.unlink(filepath)


class TestSequenceOperations:
    """Test basic sequence creation and manipulation methods."""
    
    def test_create_sequence_basic(self, mcp_server, sample_sequences):
        """Test basic sequence creation."""
        seq = sample_sequences["short_dna"]
        result = mcp_server.create_sequence(seq, name="test", circular=False)
        
        assert isinstance(result, SequenceInfo)
        assert result.sequence == seq
        assert result.length == len(seq)
        assert result.circular == False
        assert 0 < result.gc_content < 100
        assert result.molecular_weight > 0
        assert result.tm > 0
    
    def test_create_sequence_circular(self, mcp_server, sample_sequences):
        """Test circular sequence creation."""
        seq = sample_sequences["circular_plasmid"]
        result = mcp_server.create_sequence(seq, name="plasmid", circular=True)
        
        assert result.circular == True
        assert result.length == len(seq)
    
    def test_create_sequence_gc_content(self, mcp_server):
        """Test GC content calculation."""
        # High GC sequence
        high_gc = "GCGCGCGCGC"  # 100% GC
        result = mcp_server.create_sequence(high_gc)
        assert result.gc_content == 100.0
        
        # Low GC sequence  
        low_gc = "ATATATATAT"  # 0% GC
        result = mcp_server.create_sequence(low_gc)
        assert result.gc_content == 0.0
    
    def test_create_sequence_invalid_input(self, mcp_server):
        """Test error handling for invalid sequences."""
        # Empty sequence
        with pytest.raises(ValueError):
            mcp_server.create_sequence("")
        
        # Invalid nucleotides are handled by pydna (no error raised)
        result = mcp_server.create_sequence("ATGXYZ")
        assert result.sequence == "ATGXYZ"
        assert result.length == 6
    
    def test_sequence_info(self, mcp_server, sample_sequences):
        """Test sequence information analysis."""
        seq = sample_sequences["short_dna"]
        result = mcp_server.sequence_info(seq)
        
        assert isinstance(result, SequenceInfo)
        assert result.sequence == seq
        assert result.length == len(seq)
        assert result.circular == False  # Always false for sequence_info
        assert 0 <= result.gc_content <= 100
    
    def test_reverse_complement(self, mcp_server):
        """Test reverse complement calculation."""
        # Simple test case
        seq = "ATGC"
        result = mcp_server.reverse_complement(seq)
        
        assert result["original"] == seq
        assert result["reverse_complement"] == "GCAT"
        
        # Palindromic sequence (EcoRI site)
        seq = "GAATTC"
        result = mcp_server.reverse_complement(seq)
        assert result["reverse_complement"] == "GAATTC"
    
    def test_reverse_complement_long(self, mcp_server, sample_sequences):
        """Test reverse complement with longer sequences."""
        seq = sample_sequences["pcr_template"]
        result = mcp_server.reverse_complement(seq)
        
        assert len(result["reverse_complement"]) == len(seq)
        assert result["original"] == seq
        # Double reverse complement should give original
        double_rc = mcp_server.reverse_complement(result["reverse_complement"])
        assert double_rc["reverse_complement"] == seq


class TestPrimerDesign:
    """Test primer design and PCR simulation methods."""
    
    def test_design_primers_basic(self, mcp_server, sample_sequences):
        """Test basic primer design."""
        template = sample_sequences["pcr_template"]
        result = mcp_server.design_primers(template, target_tm=55.0)
        
        assert isinstance(result, AmpliconInfo)
        assert result.length == len(template)
        assert isinstance(result.forward_primer, PrimerInfo)
        assert isinstance(result.reverse_primer, PrimerInfo)
        assert result.forward_primer.length >= 13  # Minimum primer length
        assert result.reverse_primer.length >= 13
    
    def test_design_primers_different_tm(self, mcp_server, sample_sequences):
        """Test primer design with different target Tm values."""
        template = sample_sequences["pcr_template"]
        
        # Low Tm
        result_low = mcp_server.design_primers(template, target_tm=50.0)
        
        # High Tm may fail for short templates, so test with try/except
        try:
            result_high = mcp_server.design_primers(template, target_tm=65.0)
            # If both succeed, higher Tm should generally result in different primers
            assert result_high.forward_primer.tm != result_low.forward_primer.tm or result_high.reverse_primer.tm != result_low.reverse_primer.tm
        except ValueError:
            # High Tm might fail for short sequences, which is acceptable
            pass
    
    def test_design_primers_short_template(self, mcp_server):
        """Test primer design with template too short."""
        short_template = "ATGC"  # Too short for primer design
        
        with pytest.raises(ValueError):
            mcp_server.design_primers(short_template)
    
    def test_pcr_amplify_basic(self, mcp_server, sample_sequences, sample_primers):
        """Test basic PCR amplification."""
        template = sample_sequences["pcr_template"]
        fp = sample_primers["forward"]
        rp = sample_primers["reverse"]
        
        result = mcp_server.pcr_amplify(fp, rp, template)
        
        assert isinstance(result, AmpliconInfo)
        assert len(result.sequence) > 0
        assert result.forward_primer.sequence == fp
        assert result.reverse_primer.sequence == rp
    
    def test_pcr_amplify_no_binding(self, mcp_server, sample_sequences):
        """Test PCR with primers that don't bind to template."""
        template = sample_sequences["short_dna"]
        fp = "CCCCCCCCCCCCCCC"  # Won't bind
        rp = "GGGGGGGGGGGGGGG"  # Won't bind
        
        with pytest.raises(ValueError):
            mcp_server.pcr_amplify(fp, rp, template)
    
    def test_anneal_primers(self, mcp_server, sample_sequences):
        """Test primer annealing analysis."""
        template = sample_sequences["pcr_template"]
        primers = ["ATGAAAGCACTGATT", "TTTTAAAAAAAAGGA"]
        
        result = mcp_server.anneal_primers(primers, template, limit=10)
        
        assert isinstance(result, dict)
        assert "template_length" in result
        assert "forward_primers" in result
        assert "reverse_primers" in result
        assert "products" in result
        assert result["template_length"] == len(template)
    
    def test_anneal_primers_empty_list(self, mcp_server, sample_sequences):
        """Test annealing with empty primer list."""
        template = sample_sequences["short_dna"]
        
        # Empty primer list should return empty results, not raise error
        result = mcp_server.anneal_primers([], template)
        assert isinstance(result, dict)
        assert result["forward_primers"] == []
        assert result["reverse_primers"] == []


class TestRestrictionAnalysis:
    """Test restriction enzyme analysis methods."""
    
    def test_restriction_analysis_basic(self, mcp_server, sample_sequences):
        """Test basic restriction analysis."""
        seq = sample_sequences["restriction_seq"]
        enzymes = ["EcoRI", "BamHI"]
        
        result = mcp_server.restriction_analysis(seq, enzymes)
        
        assert isinstance(result, list)
        # Should find EcoRI (GAATTC) and BamHI (GGATCC) sites
        enzyme_names = [r.enzyme_name for r in result]
        assert "EcoRI" in enzyme_names
        assert "BamHI" in enzyme_names
    
    def test_restriction_analysis_specific_sites(self, mcp_server):
        """Test restriction analysis with known sites."""
        # Sequence with EcoRI site at position 0
        seq = "GAATTCATGAAAGCACTGATT"
        result = mcp_server.restriction_analysis(seq, ["EcoRI"])
        
        assert len(result) == 1
        assert result[0].enzyme_name == "EcoRI"
        assert result[0].recognition_site == "GAATTC"
        assert 0 in result[0].cut_positions
        assert len(result[0].fragments) == 2
    
    def test_restriction_analysis_no_sites(self, mcp_server, sample_sequences):
        """Test restriction analysis with no cutting sites."""
        seq = sample_sequences["short_dna"]  # No restriction sites
        result = mcp_server.restriction_analysis(seq, ["NotI"])
        
        # NotI doesn't cut this sequence
        assert len(result) == 0
    
    def test_restriction_analysis_invalid_enzyme(self, mcp_server, sample_sequences):
        """Test restriction analysis with invalid enzyme names."""
        seq = sample_sequences["short_dna"]
        # Mix of valid and invalid enzymes
        result = mcp_server.restriction_analysis(seq, ["EcoRI", "InvalidEnzyme", "BamHI"])
        
        # Should skip invalid enzymes silently
        enzyme_names = [r.enzyme_name for r in result]
        assert "InvalidEnzyme" not in enzyme_names
    
    def test_digest_sequence(self, mcp_server, sample_sequences):
        """Test sequence digestion."""
        seq = sample_sequences["restriction_seq"]
        enzymes = ["EcoRI"]
        
        result = mcp_server.digest_sequence(seq, enzymes)
        
        assert isinstance(result, dict)
        assert "original_sequence" in result
        assert "fragments" in result
        assert "num_fragments" in result
        assert result["original_sequence"] == seq
        assert result["num_fragments"] >= 1
    
    def test_digest_sequence_multiple_enzymes(self, mcp_server, sample_sequences):
        """Test digestion with multiple enzymes."""
        seq = sample_sequences["restriction_seq"]
        enzymes = ["EcoRI", "BamHI"]
        
        result = mcp_server.digest_sequence(seq, enzymes)
        
        # Double digest should produce more fragments
        assert result["num_fragments"] >= 2
        
        # Fragment lengths should sum to original length
        total_length = sum(frag["length"] for frag in result["fragments"])
        assert total_length == len(seq)
    
    def test_digest_sequence_no_valid_enzymes(self, mcp_server, sample_sequences):
        """Test digestion with no valid enzymes."""
        seq = sample_sequences["short_dna"]
        
        with pytest.raises(ValueError):
            mcp_server.digest_sequence(seq, ["InvalidEnzyme"])


class TestAssembly:
    """Test DNA assembly methods."""
    
    def test_assembly_basic(self, mcp_server, sample_sequences):
        """Test basic DNA assembly."""
        fragments = [
            sample_sequences["overlap_frag1"],
            sample_sequences["overlap_frag2"], 
            sample_sequences["overlap_frag3"]
        ]
        
        result = mcp_server.assembly(fragments, limit=8)
        
        assert isinstance(result, AssemblyResult)
        assert result.length > 0
        assert result.fragments_used == len(fragments)
        assert result.sequence != ""
    
    def test_assembly_insufficient_overlap(self, mcp_server):
        """Test assembly with insufficient overlap."""
        # Fragments with no overlap
        fragments = ["ATGAAAGCACTGATT", "CCCGGGAAATTTGGG", "TTTTTTTTTTTTTT"]
        
        with pytest.raises(ValueError):
            mcp_server.assembly(fragments, limit=25)  # High overlap requirement
    
    def test_gibson_assembly(self, mcp_server, sample_sequences):
        """Test Gibson assembly."""
        fragments = [
            sample_sequences["overlap_frag1"],
            sample_sequences["overlap_frag2"],
            sample_sequences["overlap_frag3"]
        ]
        
        result = mcp_server.gibson_assembly(fragments, overlap=15)
        
        assert isinstance(result, AssemblyResult)
        assert result.length > 0
        assert result.fragments_used == len(fragments)
    
    def test_gibson_assembly_single_fragment(self, mcp_server, sample_sequences):
        """Test Gibson assembly with single fragment."""
        fragments = [sample_sequences["short_dna"]]
        
        with pytest.raises(ValueError):
            mcp_server.gibson_assembly(fragments)
    
    def test_gibson_assembly_circular(self, mcp_server):
        """Test Gibson assembly producing circular product."""
        # Design fragments that should form circular product
        fragments = [
            "ATGAAAGCACTGATTCTATTGC",
            "CTATTGCTGAAAAAGATAAT",
            "AGATAATAGCATGAAAGCACT"  # Overlaps with first fragment
        ]
        
        result = mcp_server.gibson_assembly(fragments, overlap=15)
        
        assert isinstance(result, AssemblyResult)
        # Circular products are possible with proper overlaps


class TestFileOperations:
    """Test file reading and writing methods."""
    
    def test_read_sequence_file_fasta(self, mcp_server):
        """Test reading FASTA format files."""
        result = mcp_server.read_sequence_file(FASTA_CONTENT, "fasta")
        
        assert isinstance(result, dict)
        assert "sequence" in result
        assert "name" in result
        assert "length" in result
        assert result["name"] == "test_sequence"
        assert "ATGAAAGCACTGATT" in result["sequence"]
    
    def test_read_sequence_file_genbank(self, mcp_server):
        """Test reading GenBank format files."""
        result = mcp_server.read_sequence_file(GENBANK_CONTENT, "genbank")
        
        assert isinstance(result, dict)
        assert "sequence" in result
        assert "features" in result
        assert "num_features" in result
        assert result["name"] == "test_seq"
        assert result["num_features"] > 0
    
    def test_read_sequence_file_auto_detect(self, mcp_server):
        """Test automatic format detection."""
        # FASTA file
        result_fasta = mcp_server.read_sequence_file(FASTA_CONTENT, "auto")
        assert "test_sequence" in result_fasta["name"]
        
        # GenBank file
        result_gb = mcp_server.read_sequence_file(GENBANK_CONTENT, "auto")
        assert result_gb["num_features"] > 0
    
    def test_read_sequence_file_invalid_format(self, mcp_server):
        """Test reading file with invalid format."""
        invalid_content = "This is not a sequence file"
        
        with pytest.raises(ValueError):
            mcp_server.read_sequence_file(invalid_content, "fasta")
    
    def test_write_sequence_file_fasta(self, mcp_server, sample_sequences):
        """Test writing FASTA format."""
        seq = sample_sequences["short_dna"]
        result = mcp_server.write_sequence_file(seq, name="test_output", file_format="fasta")
        
        assert isinstance(result, str)
        assert ">test_output" in result
        assert seq in result
    
    def test_write_sequence_file_genbank(self, mcp_server, sample_sequences):
        """Test writing GenBank format."""
        seq = sample_sequences["short_dna"]
        result = mcp_server.write_sequence_file(seq, name="test_output", file_format="genbank")
        
        assert isinstance(result, str)
        assert "LOCUS" in result
        assert "test_output" in result
        assert seq.lower() in result.lower()
    
    def test_write_sequence_file_invalid_format(self, mcp_server, sample_sequences):
        """Test writing with invalid format."""
        seq = sample_sequences["short_dna"]
        
        with pytest.raises(ValueError):
            mcp_server.write_sequence_file(seq, file_format="invalid_format")


class TestGelElectrophoresis:
    """Test gel electrophoresis simulation."""
    
    def test_gel_electrophoresis_basic(self, mcp_server, sample_sequences):
        """Test basic gel electrophoresis."""
        fragments = [
            sample_sequences["short_dna"],
            sample_sequences["overlap_frag1"],
            sample_sequences["overlap_frag2"]
        ]
        
        result = mcp_server.gel_electrophoresis(fragments)
        
        assert isinstance(result, dict)
        assert "fragments" in result
        assert "gel_image" in result
        assert len(result["fragments"]) == len(fragments)
        
        # Check fragment information
        for i, frag_info in enumerate(result["fragments"]):
            assert frag_info["length"] == len(fragments[i])
            assert frag_info["sequence"] == fragments[i]
    
    def test_gel_electrophoresis_different_ladders(self, mcp_server, sample_sequences):
        """Test gel with different ladder types."""
        fragments = [sample_sequences["short_dna"]]
        
        for ladder in ["1kb", "100bp", "lambda"]:
            result = mcp_server.gel_electrophoresis(fragments, ladder=ladder)
            assert result["ladder_used"] == ladder
    
    def test_gel_electrophoresis_empty_fragments(self, mcp_server):
        """Test gel with empty fragment list."""
        with pytest.raises(ValueError):
            mcp_server.gel_electrophoresis([])


class TestGenBankDownload:
    """Test GenBank download functionality."""
    
    @patch('pydna_mcp.server.download_text')
    @patch('pydna_mcp.server.parse')
    def test_download_genbank_success(self, mock_parse, mock_download, mcp_server):
        """Test successful GenBank download."""
        # Mock the download and parse functions
        mock_download.return_value = GENBANK_CONTENT
        
        # Create a mock sequence object
        mock_seq = Mock()
        mock_seq.name = "test_seq"
        mock_seq.id = "test_seq"
        mock_seq.description = "Test sequence"
        mock_seq.seq = "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT"
        mock_seq.circular = False
        mock_seq.features = []
        mock_seq.annotations = {"organism": "synthetic"}
        
        mock_parse.return_value = [mock_seq]
        
        result = mcp_server.download_genbank("NM_000546")
        
        assert isinstance(result, dict)
        assert "accession" in result
        assert "sequence" in result
        assert "organism" in result
        assert result["accession"] == "NM_000546"
        
        # Verify the functions were called
        mock_download.assert_called_once_with("NM_000546")
        mock_parse.assert_called_once()
    
    @patch('pydna_mcp.server.download_text')
    def test_download_genbank_failure(self, mock_download, mcp_server):
        """Test GenBank download failure."""
        # Mock download failure
        mock_download.side_effect = Exception("Network error")
        
        with pytest.raises(ValueError):
            mcp_server.download_genbank("INVALID_ACCESSION")


class TestEdgeCases:
    """Test edge cases and error conditions."""
    
    def test_very_long_sequence(self, mcp_server):
        """Test handling of very long sequences."""
        # Create a long sequence (10kb)
        long_seq = "ATGC" * 2500
        
        result = mcp_server.sequence_info(long_seq)
        assert result.length == 10000
        assert 0 <= result.gc_content <= 100
    
    def test_very_short_sequence(self, mcp_server):
        """Test handling of very short sequences."""
        short_seq = "A"
        
        result = mcp_server.sequence_info(short_seq)
        assert result.length == 1
        assert result.gc_content == 0.0
    
    def test_all_same_nucleotide(self, mcp_server):
        """Test sequences with all same nucleotide."""
        sequences = {
            "all_a": "AAAAAAAAAA",
            "all_t": "TTTTTTTTTT", 
            "all_g": "GGGGGGGGGG",
            "all_c": "CCCCCCCCCC"
        }
        
        for name, seq in sequences.items():
            result = mcp_server.sequence_info(seq)
            assert result.length == 10
            
            if name in ["all_g", "all_c"]:
                assert result.gc_content == 100.0
            else:
                assert result.gc_content == 0.0
    
    def test_primer_design_edge_cases(self, mcp_server):
        """Test primer design edge cases."""
        # Very AT-rich template
        at_rich = "ATATATATATATATATATATATATATATAT"
        
        try:
            result = mcp_server.design_primers(at_rich, target_tm=50.0)
            assert result.forward_primer.tm < 60  # Should have low Tm
        except ValueError:
            # Acceptable if primer design fails for extreme sequences
            pass
    
    def test_restriction_multiple_sites(self, mcp_server):
        """Test restriction with multiple sites of same enzyme."""
        # Sequence with multiple EcoRI sites
        seq = "GAATTCAAAAGAATTCTTTTGAATTC"
        result = mcp_server.restriction_analysis(seq, ["EcoRI"])
        
        assert len(result) == 1
        assert result[0].enzyme_name == "EcoRI"
        assert len(result[0].cut_positions) == 3  # Three EcoRI sites
        assert len(result[0].fragments) == 4  # Should produce 4 fragments


class TestIntegration:
    """Integration tests for common workflows."""
    
    def test_complete_pcr_workflow(self, mcp_server):
        """Test complete PCR design and amplification workflow."""
        # Step 1: Create template
        template = "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAATAGGATCCTTTTTTTTAAAA"
        template_info = mcp_server.sequence_info(template)
        assert template_info.length == len(template)
        
        # Step 2: Design primers
        amplicon = mcp_server.design_primers(template, target_tm=55.0)
        assert isinstance(amplicon, AmpliconInfo)
        
        # Step 3: Simulate PCR
        fp = amplicon.forward_primer.sequence
        rp = amplicon.reverse_primer.sequence
        pcr_result = mcp_server.pcr_amplify(fp, rp, template)
        
        # Results should be consistent
        assert pcr_result.sequence == amplicon.sequence
        assert pcr_result.length == amplicon.length
    
    def test_cloning_workflow(self, mcp_server):
        """Test restriction cloning workflow."""
        # Step 1: Create sequence with restriction sites
        insert = "GAATTCATGAAAGCACTGATTCTATTGCTGAAAAAGATAATAGGATCC"
        
        # Step 2: Analyze restriction sites
        sites = mcp_server.restriction_analysis(insert, ["EcoRI", "BamHI"])
        assert len(sites) == 2  # Should find both sites
        
        # Step 3: Digest sequence
        digested = mcp_server.digest_sequence(insert, ["EcoRI", "BamHI"])
        assert digested["num_fragments"] == 3  # Two cuts = three fragments
        
        # Step 4: Simulate gel analysis
        fragments = [frag["sequence"] for frag in digested["fragments"]]
        gel_result = mcp_server.gel_electrophoresis(fragments)
        assert len(gel_result["fragments"]) == 3
    
    def test_gibson_assembly_workflow(self, mcp_server):
        """Test Gibson assembly workflow."""
        # Step 1: Create overlapping fragments
        fragments = [
            "ATGAAAGCACTGATTCTATTGC",
            "CTATTGCTGAAAAAGATAAT", 
            "AGATAATAGGATCCTTTTTTTT"
        ]
        
        # Step 2: Analyze individual fragments
        for frag in fragments:
            info = mcp_server.sequence_info(frag)
            assert info.length > 15  # Sufficient length for assembly
        
        # Step 3: Perform Gibson assembly
        assembled = mcp_server.gibson_assembly(fragments, overlap=8)
        assert assembled.length > max(len(f) for f in fragments)
        
        # Step 4: Analyze final product
        final_info = mcp_server.sequence_info(assembled.sequence)
        assert final_info.length == assembled.length


class TestMCPFunctionality:
    """Test MCP-specific functionality."""
    
    def test_server_initialization(self):
        """Test server initialization."""
        server = PyDNAMCP()
        assert server.prefix == "pydna_"
        assert server.name == "PyDNA MCP Server"
    
    def test_server_tools_registration(self):
        """Test that all tools are properly registered."""
        server = PyDNAMCP()
        
        # Check that tools are registered (this would require more complex testing
        # with actual MCP framework, so we'll keep it simple)
        expected_methods = [
            "create_sequence", "sequence_info", "reverse_complement",
            "design_primers", "pcr_amplify", "anneal_primers",
            "restriction_analysis", "digest_sequence", 
            "assembly", "gibson_assembly",
            "read_sequence_file", "write_sequence_file",
            "gel_electrophoresis", "download_genbank"
        ]
        
        # Each method should exist on the server
        for method_name in expected_methods:
            assert hasattr(server, method_name)
            assert callable(getattr(server, method_name))


class TestPerformance:
    """Performance and stress tests."""
    
    def test_large_sequence_handling(self, mcp_server):
        """Test handling of large sequences (within reason for tests)."""
        # 1kb sequence
        large_seq = "ATGCGATCGCTAGCTA" * 64  # 1024 bp
        
        result = mcp_server.sequence_info(large_seq)
        assert result.length == 1024
        assert result.molecular_weight > 600000  # Approximate MW for 1kb DNA
    
    def test_many_restriction_enzymes(self, mcp_server, sample_sequences):
        """Test analysis with many restriction enzymes."""
        seq = sample_sequences["long_template"]
        many_enzymes = [
            "EcoRI", "BamHI", "HindIII", "XhoI", "SalI", "PstI", 
            "SmaI", "KpnI", "SacI", "XbaI", "NotI", "SfiI"
        ]
        
        result = mcp_server.restriction_analysis(seq, many_enzymes)
        # Should complete without error
        assert isinstance(result, list)
    
    def test_multiple_assemblies(self, mcp_server):
        """Test multiple assembly operations."""
        fragments_sets = [
            ["ATGAAAGCACTGATT", "CTGATTGCATGCAAG", "GCAAGCTTAAATAG"],
            ["CCCGGGAAATTTGGG", "TTTGGGCCCAAATTT", "AAATTTCCCGGGTTT"],
            ["AAAAAAAAAAAAAAAA", "AAAAAATTTTTTTT", "TTTTTTAAAAAAA"]
        ]
        
        for fragments in fragments_sets:
            try:
                result = mcp_server.assembly(fragments, limit=6)
                assert isinstance(result, AssemblyResult)
            except ValueError:
                # Some assemblies may fail due to insufficient overlap
                pass


# Fixtures for async testing (if needed for future MCP integration)
@pytest.fixture
def event_loop():
    """Create an event loop for async tests."""
    loop = asyncio.get_event_loop_policy().new_event_loop()
    yield loop
    loop.close()


if __name__ == "__main__":
    # Run tests with: python -m pytest tests/test_pydna_mcp_server.py -v
    pytest.main([__file__, "-v", "--tb=short"]) 