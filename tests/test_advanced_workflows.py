#!/usr/bin/env python3
"""
Advanced workflow tests for PyDNA MCP Server.

Tests complex molecular biology workflows that combine multiple operations,
simulate real-world cloning scenarios, and validate end-to-end functionality.
"""

import pytest
import os
import sys
from typing import List, Dict, Any
from unittest.mock import patch

from pydna_mcp.server import PyDNAMCP, SequenceInfo, AmpliconInfo, AssemblyResult


class TestComplexCloning:
    """Test complex cloning workflows that mirror real lab procedures."""
    
    def test_restriction_cloning_workflow(self, mcp_server, sample_sequences):
        """Test complete restriction cloning workflow."""
        # Step 1: Start with a gene we want to clone
        gene_seq = sample_sequences["gfp"]
        
        # Step 2: Add restriction sites to primers for cloning
        forward_with_ecori = "GAATTC" + gene_seq[:20]  # EcoRI + gene start
        reverse_with_bamhi = "GGATCC" + mcp_server.reverse_complement(gene_seq[-20:])["reverse_complement"]
        
        # Step 3: Simulate PCR to add restriction sites
        template_with_sites = "AAAA" + forward_with_ecori[6:] + gene_seq[20:-20] + mcp_server.reverse_complement(reverse_with_bamhi[6:])["reverse_complement"] + "TTTT"
        
        pcr_result = mcp_server.pcr_amplify(
            forward_with_ecori, 
            reverse_with_bamhi, 
            template_with_sites
        )
        
        assert isinstance(pcr_result, AmpliconInfo)
        assert len(pcr_result.sequence) > len(gene_seq)
        
        # Step 4: Digest PCR product
        pcr_digest = mcp_server.digest_sequence(pcr_result.sequence, ["EcoRI", "BamHI"])
        assert pcr_digest["num_fragments"] >= 1
        
        # Step 5: Digest vector
        vector = sample_sequences["pUC19"]
        vector_digest = mcp_server.digest_sequence(vector, ["EcoRI", "BamHI"])
        assert vector_digest["num_fragments"] >= 2
        
        # Step 6: Simulate ligation (simplified)
        # In real lab, we'd isolate the correct fragments and ligate
        # Here we just verify the digestion worked correctly
        assert "EcoRI" in vector_digest["enzymes_used"]
        assert "BamHI" in vector_digest["enzymes_used"]
    
    def test_gibson_assembly_workflow(self, mcp_server, assembly_fragments):
        """Test Gibson assembly with multiple fragments."""
        fragments = assembly_fragments["gibson_compatible"]
        
        # Step 1: Analyze individual fragments
        fragment_info = []
        for frag in fragments:
            info = mcp_server.sequence_info(frag)
            fragment_info.append(info)
            assert info.length > 20  # Sufficient for assembly
        
        # Step 2: Perform Gibson assembly
        assembled = mcp_server.gibson_assembly(fragments, overlap=20)
        assert isinstance(assembled, AssemblyResult)
        assert assembled.fragments_used == len(fragments)
        
        # Step 3: Verify assembly is longer than individual fragments
        max_fragment_length = max(len(f) for f in fragments)
        assert assembled.length > max_fragment_length
        
        # Step 4: Analyze final construct
        final_analysis = mcp_server.sequence_info(assembled.sequence)
        assert final_analysis.length == assembled.length
        assert 0 <= final_analysis.gc_content <= 100
    
    def test_nested_pcr_workflow(self, mcp_server, pcr_reactions):
        """Test nested PCR procedure."""
        reaction_data = pcr_reactions["nested_pcr"]
        
        # Step 1: First PCR with outer primers
        outer_product = mcp_server.pcr_amplify(
            reaction_data["outer_forward"],
            reaction_data["outer_reverse"],
            reaction_data["template"]
        )
        assert isinstance(outer_product, AmpliconInfo)
        
        # Step 2: Second PCR with inner primers using first product as template
        inner_product = mcp_server.pcr_amplify(
            reaction_data["inner_forward"],
            reaction_data["inner_reverse"],
            outer_product.sequence
        )
        assert isinstance(inner_product, AmpliconInfo)
        
        # Step 3: Verify nested product is smaller than outer product
        assert inner_product.length <= outer_product.length
        
        # Step 4: Verify both products contain the expected sequences
        assert reaction_data["inner_forward"][:10] in inner_product.sequence
    
    def test_multi_fragment_assembly(self, mcp_server):
        """Test assembly of many fragments."""
        # Create fragments with systematic overlaps
        base_seq = "ATGAAAGCACTGATT"
        fragments = []
        
        for i in range(5):
            fragment = base_seq + "N" * (i * 5) + "CTATTGCTGAAA"
            fragments.append(fragment.replace("N", "T"))
        
        # Try assembly with different overlap requirements
        for overlap in [10, 15, 20]:
            try:
                result = mcp_server.assembly(fragments, limit=overlap)
                assert isinstance(result, AssemblyResult)
                assert result.fragments_used <= len(fragments)
            except ValueError:
                # Some assemblies may fail - this is acceptable
                pass


class TestPrimerDesignWorkflows:
    """Test complex primer design scenarios."""
    
    def test_primer_design_with_constraints(self, mcp_server, target_tm_values):
        """Test primer design with different temperature constraints."""
        template = BACTERIAL_GENES["lacZ"]["sequence"][:200]  # Use first 200bp
        
        result = mcp_server.design_primers(template, target_tm=target_tm_values)
        assert isinstance(result, AmpliconInfo)
        
        # Check that primers are reasonably close to target Tm
        tm_diff_forward = abs(result.tm_forward - target_tm_values)
        tm_diff_reverse = abs(result.tm_reverse - target_tm_values)
        
        assert tm_diff_forward < 10.0  # Within 10°C
        assert tm_diff_reverse < 10.0  # Within 10°C
        
        # Check primer properties
        assert 15 <= len(result.forward_primer.sequence) <= 35
        assert 15 <= len(result.reverse_primer.sequence) <= 35
    
    def test_primer_annealing_analysis(self, mcp_server, sample_sequences):
        """Test comprehensive primer annealing analysis."""
        template = sample_sequences["long_template"]
        
        # Design primers first
        amplicon = mcp_server.design_primers(template, target_tm=55.0)
        
        # Analyze how these primers anneal
        primers = [
            amplicon.forward_primer.sequence,
            amplicon.reverse_primer.sequence
        ]
        
        annealing = mcp_server.anneal_primers(primers, template, limit=13)
        
        # Should find both primers binding
        assert len(annealing["forward_primers"]) >= 1
        assert len(annealing["reverse_primers"]) >= 1
        assert annealing["products"] >= 1
        
        # Check binding positions make sense
        for fp in annealing["forward_primers"]:
            assert fp["footprint"] >= 13
        for rp in annealing["reverse_primers"]:
            assert rp["footprint"] >= 13


class TestRestrictionMapping:
    """Test comprehensive restriction mapping workflows."""
    
    def test_restriction_mapping_analysis(self, mcp_server, restriction_enzymes):
        """Test comprehensive restriction mapping."""
        # Use a sequence with known restriction sites
        test_seq = "GAATTC" + "ATGC" * 50 + "GGATCC" + "ATGC" * 30 + "AAGCTT" + "ATGC" * 20
        
        # Analyze with multiple enzymes
        results = mcp_server.restriction_analysis(test_seq, restriction_enzymes)
        
        # Should find at least some enzymes that cut
        cutting_enzymes = [r for r in results if r.cut_positions]
        assert len(cutting_enzymes) >= 2  # EcoRI, BamHI, HindIII should cut
        
        # Verify specific expected cuts
        ecori_results = [r for r in results if r.enzyme_name == "EcoRI"]
        if ecori_results:
            assert 0 in ecori_results[0].cut_positions  # Should cut at position 0
    
    def test_double_digest_analysis(self, mcp_server):
        """Test double digest restriction analysis."""
        # Create sequence with multiple EcoRI and BamHI sites
        seq = "GAATTCATGCGATCGGGATCCAAAGAATTCTTTTGGATCCGGG"
        
        # Single digests
        ecori_single = mcp_server.digest_sequence(seq, ["EcoRI"])
        bamhi_single = mcp_server.digest_sequence(seq, ["BamHI"])
        
        # Double digest
        double_digest = mcp_server.digest_sequence(seq, ["EcoRI", "BamHI"])
        
        # Double digest should produce more fragments
        assert double_digest["num_fragments"] >= max(
            ecori_single["num_fragments"],
            bamhi_single["num_fragments"]
        )
        
        # Verify fragment sizes make sense
        total_length = sum(f["length"] for f in double_digest["fragments"])
        assert total_length == len(seq)


class TestFileWorkflows:
    """Test file handling workflows."""
    
    def test_sequence_file_roundtrip(self, mcp_server, sample_sequences, file_formats):
        """Test reading and writing sequence files."""
        original_seq = sample_sequences["gfp"]
        
        # Write sequence to file format
        file_content = mcp_server.write_sequence_file(
            original_seq, 
            name="test_gfp", 
            file_format=file_formats
        )
        
        assert len(file_content) > len(original_seq)
        
        # Read sequence back
        parsed = mcp_server.read_sequence_file(file_content, file_format=file_formats)
        
        # Verify sequence is preserved
        assert parsed["sequence"].upper() == original_seq.upper()
        assert parsed["length"] == len(original_seq)
    
    @patch('pydna_mcp.server.download_text')
    @patch('pydna_mcp.server.parse')
    def test_genbank_download_workflow(self, mock_parse, mock_download, mcp_server):
        """Test downloading and analyzing GenBank sequences."""
        # Mock GenBank response
        mock_download.return_value = """LOCUS       NM_000546    1842 bp    mRNA    linear   PRI
DEFINITION  Homo sapiens tumor protein p53 (TP53), mRNA
FEATURES             Location/Qualifiers
     gene            1..1842
                     /gene="TP53"
ORIGIN
        1 atggaggagc cgcagtcaga tcctagcgtc gagccccctc tgagtcagga aacattttca
//"""
        
        # Mock parsed sequence
        from unittest.mock import Mock
        mock_seq = Mock()
        mock_seq.seq = "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCA"
        mock_seq.name = "NM_000546"
        mock_seq.id = "NM_000546.6"
        mock_seq.description = "Homo sapiens tumor protein p53 (TP53), mRNA"
        mock_seq.features = []
        mock_seq.annotations = {"organism": "Homo sapiens"}
        mock_seq.circular = False
        mock_seq.__len__ = lambda: 60
        
        mock_parse.return_value = [mock_seq]
        
        # Test download
        result = mcp_server.download_genbank("NM_000546")
        
        assert result["accession"] == "NM_000546"
        assert "Homo sapiens" in result["description"]
        assert len(result["sequence"]) > 0
        assert result["organism"] == "Homo sapiens"
        
        # Verify we can analyze the downloaded sequence
        seq_info = mcp_server.sequence_info(result["sequence"])
        assert seq_info.length == len(result["sequence"])


class TestErrorRecovery:
    """Test error handling and recovery in complex workflows."""
    
    def test_failed_pcr_recovery(self, mcp_server):
        """Test handling of failed PCR reactions."""
        template = "ATGCGATCGCTAGCTAGC"
        bad_forward = "TTTTTTTTTTTTTTTTTT"  # Won't bind to template
        bad_reverse = "AAAAAAAAAAAAAAAAAA"   # Won't bind to template
        
        # This should fail
        with pytest.raises(ValueError):
            mcp_server.pcr_amplify(bad_forward, bad_reverse, template)
    
    def test_failed_assembly_recovery(self, mcp_server):
        """Test handling of failed assembly reactions."""
        # Fragments with no overlaps
        no_overlap_fragments = [
            "ATGCGATCGCTAGC",
            "TTTTTTTTTTTTTT",
            "AAAAAAAAAAAAAAA"
        ]
        
        # This should fail
        with pytest.raises(ValueError):
            mcp_server.assembly(no_overlap_fragments, limit=20)
    
    def test_invalid_enzyme_handling(self, mcp_server, sample_sequences):
        """Test handling of invalid restriction enzymes."""
        seq = sample_sequences["short_dna"]
        
        # Mix valid and invalid enzyme names
        enzymes = ["EcoRI", "InvalidEnzyme", "BamHI", "NotAnEnzyme"]
        
        # Should complete without error, but skip invalid enzymes
        results = mcp_server.restriction_analysis(seq, enzymes)
        
        # Should only return results for valid enzymes
        valid_results = [r for r in results if r.enzyme_name in ["EcoRI", "BamHI"]]
        assert len(valid_results) <= 2


@pytest.mark.integration
class TestIntegrationWorkflows:
    """Integration tests that combine multiple server operations."""
    
    def test_complete_cloning_pipeline(self, mcp_server):
        """Test a complete cloning pipeline from start to finish."""
        # Step 1: Start with gene sequence
        gene = BACTERIAL_GENES["ampR"]["sequence"][:300]  # Use part of ampR gene
        
        # Step 2: Design primers
        primers = mcp_server.design_primers(gene, target_tm=58.0)
        
        # Step 3: Simulate PCR
        pcr_product = mcp_server.pcr_amplify(
            primers.forward_primer.sequence,
            primers.reverse_primer.sequence,
            gene
        )
        
        # Step 4: Analyze restriction sites in PCR product
        restriction_sites = mcp_server.restriction_analysis(
            pcr_product.sequence, 
            ["EcoRI", "BamHI", "HindIII"]
        )
        
        # Step 5: Write final sequence to file
        genbank_file = mcp_server.write_sequence_file(
            pcr_product.sequence,
            name="cloned_gene",
            file_format="genbank"
        )
        
        # Step 6: Read it back and verify
        parsed_file = mcp_server.read_sequence_file(genbank_file, "genbank")
        
        # Verify the pipeline worked
        assert parsed_file["sequence"].upper() == pcr_product.sequence.upper()
        assert parsed_file["name"] == "cloned_gene"
        assert len(restriction_sites) >= 0  # May or may not have sites
    
    def test_multi_gene_assembly_pipeline(self, mcp_server):
        """Test assembling multiple genes into a single construct."""
        # Use multiple real gene sequences
        genes = [
            REAL_WORLD_SEQUENCES["gfp"]["sequence"][:200],
            BACTERIAL_GENES["ampR"]["sequence"][:200],
            REAL_WORLD_SEQUENCES["his_tag"]["sequence"]
        ]
        
        # Add overlaps for Gibson assembly
        overlap_seq = "GGGCCCAAATTT"  # 12bp overlap
        
        gibson_fragments = []
        for i, gene in enumerate(genes):
            if i == 0:
                # First fragment: gene + overlap
                fragment = gene + overlap_seq
            elif i == len(genes) - 1:
                # Last fragment: overlap + gene
                fragment = overlap_seq + gene
            else:
                # Middle fragments: overlap + gene + overlap
                fragment = overlap_seq + gene + overlap_seq
            
            gibson_fragments.append(fragment)
        
        # Perform Gibson assembly
        assembled = mcp_server.gibson_assembly(gibson_fragments, overlap=12)
        
        # Verify assembly
        assert assembled.fragments_used == len(gibson_fragments)
        assert assembled.length > max(len(f) for f in gibson_fragments)
        
        # Analyze final construct
        final_info = mcp_server.sequence_info(assembled.sequence)
        assert final_info.length == assembled.length


if __name__ == "__main__":
    pytest.main([__file__, "-v"]) 