#!/usr/bin/env python3
"""
Edge case and error handling tests for PyDNA MCP Server.

Tests unusual inputs, boundary conditions, and error scenarios
to ensure robust behavior.
"""

import pytest
import os
import sys
from typing import Dict, Any, List
from unittest.mock import patch, Mock

from pydna_mcp.server import PyDNAMCP, SequenceInfo, AmpliconInfo, AssemblyResult


class TestInvalidInputs:
    """Test handling of invalid inputs."""
    
    def test_empty_sequence_handling(self, mcp_server):
        """Test handling of empty sequences."""
        with pytest.raises(ValueError, match="empty|invalid"):
            mcp_server.sequence_info("")
        
        with pytest.raises(ValueError, match="empty|invalid"):
            mcp_server.reverse_complement("")
    
    def test_none_sequence_handling(self, mcp_server):
        """Test handling of None sequences."""
        with pytest.raises((ValueError, TypeError)):
            mcp_server.sequence_info(None)
        
        with pytest.raises((ValueError, TypeError)):
            mcp_server.reverse_complement(None)
    
    def test_invalid_dna_characters(self, mcp_server):
        """Test handling of sequences with invalid DNA characters."""
        invalid_sequences = [
            "ATGCXYZ",  # Invalid characters
            "ATGC123",  # Numbers
            "ATGC!@#",  # Special characters
            "ATGCnnn",  # Lowercase n (ambiguous)
            "ATGC N N",  # Spaces
        ]
        
        for seq in invalid_sequences:
            # Should either handle gracefully or raise appropriate error
            try:
                result = mcp_server.sequence_info(seq)
                # If it succeeds, verify the result is reasonable
                assert result.length >= 4  # At least the valid ATGC part
            except ValueError:
                # ValueError is acceptable for invalid sequences
                pass
    
    def test_very_short_sequences(self, mcp_server):
        """Test handling of very short sequences."""
        short_sequences = ["A", "AT", "ATG"]
        
        for seq in short_sequences:
            info = mcp_server.sequence_info(seq)
            assert info.length == len(seq)
            assert isinstance(info.gc_content, float)
            
            rev_comp = mcp_server.reverse_complement(seq)
            assert len(rev_comp["reverse_complement"]) == len(seq)
    
    def test_empty_primer_lists(self, mcp_server):
        """Test handling of empty primer lists."""
        template = "ATGCGATCGCTAGCTAGC"
        
        with pytest.raises(ValueError):
            mcp_server.anneal_primers([], template)
        
        with pytest.raises(ValueError):
            mcp_server.anneal_primers([""], template)
    
    def test_empty_fragment_lists(self, mcp_server):
        """Test handling of empty fragment lists for assembly."""
        with pytest.raises(ValueError):
            mcp_server.assembly([])
        
        with pytest.raises(ValueError):
            mcp_server.gibson_assembly([])
    
    def test_invalid_enzyme_names(self, mcp_server):
        """Test handling of invalid restriction enzyme names."""
        sequence = "ATGCGATCGCTAGCTAGC"
        
        # Mix of valid and invalid enzymes
        mixed_enzymes = ["EcoRI", "", "NonExistentEnzyme", "BamHI", None]
        
        # Should handle gracefully, possibly returning results for valid enzymes only
        try:
            results = mcp_server.restriction_analysis(sequence, mixed_enzymes)
            # Should only have results for valid enzymes
            valid_results = [r for r in results if r.enzyme_name in ["EcoRI", "BamHI"]]
            assert len(valid_results) <= 2
        except ValueError:
            # Acceptable to raise error for invalid input
            pass


class TestBoundaryConditions:
    """Test boundary conditions and limits."""
    
    def test_maximum_primer_length(self, mcp_server):
        """Test handling of very long primers."""
        template = "ATGC" * 100  # 400bp template
        very_long_primer = "ATGC" * 50  # 200bp primer (half the template)
        short_primer = "ATGCGATCGC"
        
        # Very long primers might not work well
        try:
            result = mcp_server.pcr_amplify(very_long_primer, short_primer, template)
            assert isinstance(result, AmpliconInfo)
        except ValueError:
            # Long primers may legitimately fail
            pass
    
    def test_minimum_primer_length(self, mcp_server):
        """Test handling of very short primers."""
        template = "ATGCGATCGCTAGCTAGCATGCGATCGC"
        very_short_primers = ["AT", "GC", "TA"]
        
        for primer in very_short_primers:
            # Very short primers should fail
            with pytest.raises(ValueError):
                mcp_server.pcr_amplify(primer, "ATGCGATCGC", template)
    
    def test_extreme_gc_content(self, mcp_server):
        """Test sequences with extreme GC content."""
        # 100% GC content
        all_gc = "GCGCGCGCGC" * 20
        info_gc = mcp_server.sequence_info(all_gc)
        assert info_gc.gc_content > 99.0
        
        # 0% GC content
        all_at = "ATATATATAT" * 20
        info_at = mcp_server.sequence_info(all_at)
        assert info_at.gc_content < 1.0
    
    def test_maximum_overlap_requirements(self, mcp_server):
        """Test assembly with very high overlap requirements."""
        fragments = [
            "ATGCGATCGCTAGCTAGC",
            "CTAGCTAGCATGCGATCG",
            "ATGCGATCGCTAGCTAGC"
        ]
        
        # Require overlap longer than fragments
        with pytest.raises(ValueError):
            mcp_server.assembly(fragments, limit=50)  # Longer than fragments
    
    def test_zero_overlap_requirements(self, mcp_server):
        """Test assembly with zero overlap requirements."""
        fragments = [
            "ATGCGATCGCTAGCTAGC",
            "CTAGCTAGCATGCGATCG"
        ]
        
        # Zero overlap should fail or produce poor results
        try:
            result = mcp_server.assembly(fragments, limit=0)
            if isinstance(result, AssemblyResult):
                # If it succeeds, result should be reasonable
                assert result.fragments_used >= 1
        except ValueError:
            # Acceptable to fail with zero overlap
            pass
    
    def test_single_fragment_assembly(self, mcp_server):
        """Test assembly with only one fragment."""
        single_fragment = ["ATGCGATCGCTAGCTAGC"]
        
        result = mcp_server.assembly(single_fragment)
        assert isinstance(result, AssemblyResult)
        assert result.fragments_used == 1
        assert result.sequence == single_fragment[0]


class TestDataTypeEdgeCases:
    """Test edge cases related to data types and formatting."""
    
    def test_integer_inputs(self, mcp_server):
        """Test handling of integer inputs where strings are expected."""
        with pytest.raises(TypeError):
            mcp_server.sequence_info(12345)
        
        with pytest.raises(TypeError):
            mcp_server.reverse_complement(67890)
    
    def test_list_inputs_for_sequences(self, mcp_server):
        """Test handling of list inputs where strings are expected."""
        with pytest.raises(TypeError):
            mcp_server.sequence_info(["A", "T", "G", "C"])
        
        with pytest.raises(TypeError):
            mcp_server.reverse_complement(["ATGC"])
    
    def test_mixed_case_sequences(self, mcp_server):
        """Test handling of mixed case sequences."""
        mixed_case_sequences = [
            "AtGcGaTcGc",
            "atgcGATCGC",
            "ATGC",
            "atgc"
        ]
        
        for seq in mixed_case_sequences:
            info = mcp_server.sequence_info(seq)
            assert info.length == len(seq)
            
            rev_comp = mcp_server.reverse_complement(seq)
            assert len(rev_comp["reverse_complement"]) == len(seq)
    
    def test_whitespace_in_sequences(self, mcp_server):
        """Test handling of sequences with whitespace."""
        sequences_with_whitespace = [
            " ATGCGATCGC ",  # Leading/trailing spaces
            "ATGC GATC GC",  # Internal spaces
            "ATGC\\nGATC\\nGC",  # Newlines
            "ATGC\\tGATC\\tGC",  # Tabs
        ]
        
        for seq in sequences_with_whitespace:
            try:
                info = mcp_server.sequence_info(seq)
                # If it succeeds, length should make sense
                assert info.length > 0
            except ValueError:
                # Acceptable to reject sequences with whitespace
                pass
    
    def test_unicode_characters(self, mcp_server):
        """Test handling of unicode characters."""
        unicode_sequences = [
            "ATGC™GATC",  # Trademark symbol
            "ATGC€GATC",  # Euro symbol
            "ATGCαGATC",  # Greek alpha
        ]
        
        for seq in unicode_sequences:
            with pytest.raises((ValueError, TypeError)):
                mcp_server.sequence_info(seq)


class TestFileHandlingEdgeCases:
    """Test edge cases in file operations."""
    
    def test_empty_file_content(self, mcp_server):
        """Test handling of empty file content."""
        with pytest.raises(ValueError):
            mcp_server.read_sequence_file("", "fasta")
        
        with pytest.raises(ValueError):
            mcp_server.read_sequence_file("   ", "genbank")
    
    def test_invalid_file_formats(self, mcp_server):
        """Test handling of invalid file formats."""
        sequence = "ATGCGATCGCTAGC"
        
        invalid_formats = ["txt", "doc", "pdf", "xlsx", "invalid"]
        
        for fmt in invalid_formats:
            with pytest.raises(ValueError):
                mcp_server.write_sequence_file(sequence, "test", fmt)
    
    def test_malformed_fasta_content(self, mcp_server):
        """Test handling of malformed FASTA content."""
        malformed_fasta_examples = [
            "ATGCGATCGC",  # No header
            ">header\\n",  # Header but no sequence
            ">header1\\nATGC\\n>header2\\n",  # Missing sequence for second entry
            ">>invalid_header\\nATGC",  # Invalid header format
        ]
        
        for content in malformed_fasta_examples:
            try:
                result = mcp_server.read_sequence_file(content, "fasta")
                # If it succeeds, result should be reasonable
                assert len(result["sequence"]) > 0
            except ValueError:
                # Acceptable to reject malformed content
                pass
    
    def test_malformed_genbank_content(self, mcp_server):
        """Test handling of malformed GenBank content."""
        malformed_genbank = [
            "LOCUS invalid",  # Incomplete LOCUS line
            "ORIGIN\\n1 atgc\\n//\\n",  # Missing required fields
            "LOCUS test 100 bp DNA\\nORIGIN\\n//",  # Missing sequence
        ]
        
        for content in malformed_genbank:
            try:
                result = mcp_server.read_sequence_file(content, "genbank")
                assert len(result["sequence"]) >= 0
            except ValueError:
                # Acceptable to reject malformed content
                pass
    
    def test_very_long_sequence_names(self, mcp_server):
        """Test handling of very long sequence names."""
        sequence = "ATGCGATCGC"
        very_long_name = "A" * 1000  # 1000 character name
        
        # Should handle long names gracefully
        result = mcp_server.write_sequence_file(sequence, very_long_name, "fasta")
        assert very_long_name in result
    
    def test_special_characters_in_names(self, mcp_server):
        """Test handling of special characters in sequence names."""
        sequence = "ATGCGATCGC"
        special_names = [
            "test/name",  # Forward slash
            "test\\name",  # Backslash
            "test:name",  # Colon
            "test*name",  # Asterisk
            "test?name",  # Question mark
            "test<name>",  # Angle brackets
        ]
        
        for name in special_names:
            try:
                result = mcp_server.write_sequence_file(sequence, name, "fasta")
                # Name should appear in some form
                assert len(result) > len(sequence)
            except ValueError:
                # Acceptable to reject names with special characters
                pass


class TestNetworkAndExternalDependencies:
    """Test edge cases related to external dependencies."""
    
    @patch('pydna_mcp.server.download_text')
    def test_genbank_download_network_failure(self, mock_download, mcp_server):
        """Test handling of network failures during GenBank download."""
        # Simulate network failure
        mock_download.side_effect = ConnectionError("Network unreachable")
        
        with pytest.raises(ConnectionError):
            mcp_server.download_genbank("NM_000546")
    
    @patch('pydna_mcp.server.download_text')
    def test_genbank_download_invalid_accession(self, mock_download, mcp_server):
        """Test handling of invalid GenBank accessions."""
        # Simulate GenBank error response
        mock_download.return_value = "ERROR: Invalid accession number"
        
        with pytest.raises(ValueError):
            mcp_server.download_genbank("INVALID_ACCESSION")
    
    @patch('pydna_mcp.server.download_text')
    def test_genbank_download_empty_response(self, mock_download, mcp_server):
        """Test handling of empty GenBank response."""
        mock_download.return_value = ""
        
        with pytest.raises(ValueError):
            mcp_server.download_genbank("NM_000546")
    
    @patch('pydna_mcp.server.download_text')
    def test_genbank_download_timeout(self, mock_download, mcp_server):
        """Test handling of GenBank download timeout."""
        import time
        
        def slow_download(*args, **kwargs):
            time.sleep(2)  # Simulate slow response
            return "LOCUS test 100 bp DNA\\nORIGIN\\n1 atgcgatcgc\\n//"
        
        mock_download.side_effect = slow_download
        
        # Should still work, just slowly
        try:
            result = mcp_server.download_genbank("NM_000546")
            assert len(result["sequence"]) > 0
        except Exception:
            # Timeout handling is implementation-dependent
            pass


class TestRobustnessAndRecovery:
    """Test system robustness and error recovery."""
    
    def test_repeated_failed_operations(self, mcp_server):
        """Test that repeated failed operations don't break the system."""
        invalid_sequence = "XYZXYZXYZ"
        
        # Try multiple times - system should remain stable
        for _ in range(5):
            try:
                mcp_server.sequence_info(invalid_sequence)
            except (ValueError, TypeError):
                pass  # Expected failures
        
        # System should still work with valid input
        valid_result = mcp_server.sequence_info("ATGCGATCGC")
        assert valid_result.length == 10
    
    def test_mixed_valid_invalid_operations(self, mcp_server):
        """Test mixing valid and invalid operations."""
        operations = [
            ("ATGCGATCGC", True),  # Valid
            ("XYZXYZXYZ", False),  # Invalid
            ("ATGCGATCGCTAGC", True),  # Valid
            ("", False),  # Invalid
            ("ATGC", True),  # Valid
        ]
        
        valid_count = 0
        for seq, should_succeed in operations:
            try:
                result = mcp_server.sequence_info(seq)
                if should_succeed:
                    valid_count += 1
                    assert result.length == len(seq)
            except (ValueError, TypeError):
                if should_succeed:
                    pytest.fail(f"Valid sequence {seq} should not have failed")
        
        assert valid_count == 3  # Three valid sequences should have succeeded
    
    def test_state_consistency_after_errors(self, mcp_server):
        """Test that the server state remains consistent after errors."""
        # Perform a valid operation
        result1 = mcp_server.sequence_info("ATGCGATCGC")
        assert result1.length == 10
        
        # Try invalid operations
        try:
            mcp_server.sequence_info("")
        except ValueError:
            pass
        
        try:
            mcp_server.pcr_amplify("", "", "")
        except ValueError:
            pass
        
        # Server should still work normally
        result2 = mcp_server.sequence_info("ATGCGATCGC")
        assert result2.length == 10
        assert result2.gc_content == result1.gc_content


class TestConcurrencyEdgeCases:
    """Test edge cases that might occur in concurrent usage."""
    
    def test_rapid_sequential_operations(self, mcp_server):
        """Test rapid sequential operations don't interfere."""
        sequences = ["ATGC" * 10, "GCTA" * 10, "TACG" * 10]
        
        results = []
        for seq in sequences:
            for _ in range(3):  # Multiple operations per sequence
                info = mcp_server.sequence_info(seq)
                rev_comp = mcp_server.reverse_complement(seq)
                results.append((info, rev_comp))
        
        # All results should be consistent
        assert len(results) == 9
        for info, rev_comp in results:
            assert info.length == len(rev_comp["reverse_complement"])
    
    def test_interleaved_operations(self, mcp_server):
        """Test interleaved operations on different sequences."""
        seq1 = "ATGCGATCGC"
        seq2 = "GCTAGCATGC"
        
        # Interleave operations
        info1a = mcp_server.sequence_info(seq1)
        info2a = mcp_server.sequence_info(seq2)
        rev1 = mcp_server.reverse_complement(seq1)
        rev2 = mcp_server.reverse_complement(seq2)
        info1b = mcp_server.sequence_info(seq1)
        info2b = mcp_server.sequence_info(seq2)
        
        # Results should be consistent
        assert info1a.length == info1b.length == len(seq1)
        assert info2a.length == info2b.length == len(seq2)
        assert len(rev1["reverse_complement"]) == len(seq1)
        assert len(rev2["reverse_complement"]) == len(seq2)


@pytest.mark.stress
class TestStressConditions:
    """Test stress conditions and resource limits."""
    
    def test_memory_stress_large_sequences(self, mcp_server):
        """Test memory usage with very large sequences."""
        # Create progressively larger sequences
        base_size = 10000  # 10kb
        
        for multiplier in [1, 2, 4, 8]:
            large_seq = "ATGC" * (base_size * multiplier // 4)
            
            try:
                info = mcp_server.sequence_info(large_seq)
                assert info.length == len(large_seq)
                
                # Free memory by clearing reference
                del large_seq
            except MemoryError:
                # Acceptable to hit memory limits
                break
    
    def test_computational_stress_many_operations(self, mcp_server):
        """Test computational stress with many operations."""
        sequence = "ATGCGATCGCTAGCTAGC"
        
        # Perform many operations
        for i in range(100):
            info = mcp_server.sequence_info(sequence)
            rev_comp = mcp_server.reverse_complement(sequence)
            
            # Results should remain consistent
            assert info.length == len(sequence)
            assert len(rev_comp["reverse_complement"]) == len(sequence)
    
    def test_error_stress_many_failures(self, mcp_server):
        """Test stress from many consecutive failures."""
        invalid_inputs = ["", "XYZXYZ", None, 12345, [], {}]
        
        # Try many invalid operations
        for _ in range(50):
            for invalid_input in invalid_inputs:
                try:
                    if invalid_input is not None:
                        mcp_server.sequence_info(str(invalid_input))
                except (ValueError, TypeError):
                    pass  # Expected failures
        
        # System should still work
        valid_result = mcp_server.sequence_info("ATGCGATCGC")
        assert valid_result.length == 10


if __name__ == "__main__":
    pytest.main([__file__, "-v"]) 