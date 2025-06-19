#!/usr/bin/env python3
"""
Performance tests for PyDNA MCP Server.

Tests to ensure operations complete within reasonable time limits
and handle large sequences efficiently.
"""

import pytest
import time
import os
import sys
from typing import Dict, Any
from unittest.mock import patch

from pydna_mcp.server import PyDNAMCP, SequenceInfo, AmpliconInfo, AssemblyResult


class TestPerformanceBasics:
    """Test basic performance requirements."""
    
    def test_sequence_creation_performance(self, mcp_server):
        """Test sequence creation with large sequences."""
        large_sequence = "ATGC" * 2500  # 10kb sequence
        
        start_time = time.time()
        info = mcp_server.sequence_info(large_sequence)
        end_time = time.time()
        
        # Should complete within 5 seconds
        assert end_time - start_time < 5.0
        assert info.length == len(large_sequence)
        assert isinstance(info.gc_content, float)
    
    def test_reverse_complement_performance(self, mcp_server):
        """Test reverse complement calculation performance."""
        large_sequence = "ATGCGATCGCTAGCTAGC" * 1000  # ~18kb
        
        start_time = time.time()
        result = mcp_server.reverse_complement(large_sequence)
        end_time = time.time()
        
        # Should complete within 2 seconds
        assert end_time - start_time < 2.0
        assert len(result["reverse_complement"]) == len(large_sequence)
    
    def test_primer_design_performance(self, mcp_server):
        """Test primer design performance with long templates."""
        long_template = BACTERIAL_GENES["lacZ"]["sequence"] * 5  # ~5kb template
        
        start_time = time.time()
        result = mcp_server.design_primers(long_template, target_tm=55.0)
        end_time = time.time()
        
        # Should complete within 10 seconds
        assert end_time - start_time < 10.0
        assert isinstance(result, AmpliconInfo)
    
    def test_restriction_analysis_performance(self, mcp_server):
        """Test restriction analysis performance."""
        large_sequence = ("ATGCGATCGCTAGCTAGC" + 
                         "GAATTC" +  # EcoRI site
                         "ATGCGATCGCTAGCTAGC" + 
                         "GGATCC" +  # BamHI site
                         "ATGCGATCGCTAGCTAGC") * 100  # ~8kb with known sites
        
        enzymes = ["EcoRI", "BamHI", "HindIII", "XbaI", "SacI", "KpnI", "SmaI", "PstI"]
        
        start_time = time.time()
        results = mcp_server.restriction_analysis(large_sequence, enzymes)
        end_time = time.time()
        
        # Should complete within 15 seconds
        assert end_time - start_time < 15.0
        assert len(results) >= len([e for e in enzymes if e in ["EcoRI", "BamHI"]])
    
    def test_pcr_amplify_performance(self, mcp_server):
        """Test PCR amplification performance."""
        template = BACTERIAL_GENES["lacZ"]["sequence"]
        
        # Design primers first
        primers = mcp_server.design_primers(template[:500], target_tm=55.0)
        
        start_time = time.time()
        result = mcp_server.pcr_amplify(
            primers.forward_primer.sequence,
            primers.reverse_primer.sequence,
            template
        )
        end_time = time.time()
        
        # Should complete within 5 seconds
        assert end_time - start_time < 5.0
        assert isinstance(result, AmpliconInfo)


class TestScalabilityTests:
    """Test scalability with increasing data sizes."""
    
    @pytest.mark.parametrize("sequence_multiplier", [1, 2, 5, 10])
    def test_sequence_info_scaling(self, mcp_server, sequence_multiplier):
        """Test sequence info scaling with different sizes."""
        base_sequence = BACTERIAL_GENES["ampR"]["sequence"]
        test_sequence = base_sequence * sequence_multiplier
        
        start_time = time.time()
        info = mcp_server.sequence_info(test_sequence)
        end_time = time.time()
        
        # Time should scale roughly linearly, allowing some overhead
        max_expected_time = 0.5 + (sequence_multiplier * 0.1)
        assert end_time - start_time < max_expected_time
        assert info.length == len(test_sequence)
    
    @pytest.mark.parametrize("num_enzymes", [1, 5, 10, 20])
    def test_restriction_enzyme_scaling(self, mcp_server, num_enzymes):
        """Test restriction analysis scaling with enzyme count."""
        sequence = PERFORMANCE_TEST_SEQUENCES["medium_with_sites"]
        common_enzymes = ["EcoRI", "BamHI", "HindIII", "XbaI", "SacI", "KpnI", 
                         "SmaI", "PstI", "NcoI", "NotI", "XhoI", "SpeI", 
                         "ApaI", "SalI", "ClaI", "BglII", "AvrII", "NheI", 
                         "PmeI", "ScaI"]
        
        enzymes = common_enzymes[:num_enzymes]
        
        start_time = time.time()
        results = mcp_server.restriction_analysis(sequence, enzymes)
        end_time = time.time()
        
        # Time should scale with number of enzymes
        max_expected_time = 1.0 + (num_enzymes * 0.2)
        assert end_time - start_time < max_expected_time
        assert len(results) <= num_enzymes
    
    @pytest.mark.parametrize("fragment_count", [2, 3, 5, 8])
    def test_assembly_fragment_scaling(self, mcp_server, fragment_count):
        """Test assembly scaling with fragment count."""
        base_fragment = "ATGAAAGCACTGATTCTATTGCTG"
        overlap = "GGGCCCAAATTT"
        
        fragments = []
        for i in range(fragment_count):
            if i == 0:
                fragment = base_fragment + overlap
            elif i == fragment_count - 1:
                fragment = overlap + base_fragment
            else:
                fragment = overlap + base_fragment + overlap
            fragments.append(fragment)
        
        start_time = time.time()
        try:
            result = mcp_server.assembly(fragments, limit=10)
            end_time = time.time()
            
            # Time should scale with fragment count
            max_expected_time = 2.0 + (fragment_count * 0.5)
            assert end_time - start_time < max_expected_time
            assert isinstance(result, AssemblyResult)
        except ValueError:
            # Assembly might fail - that's ok for performance testing
            end_time = time.time()
            assert end_time - start_time < 5.0  # Should fail quickly


class TestMemoryEfficiency:
    """Test memory usage patterns."""
    
    def test_large_sequence_memory_usage(self, mcp_server):
        """Test that large sequences don't cause memory issues."""
        # Create sequences of increasing size
        sizes = [1000, 5000, 10000, 50000]  # Up to 50kb
        
        for size in sizes:
            large_seq = "ATGC" * (size // 4)
            
            # These operations should complete without memory errors
            info = mcp_server.sequence_info(large_seq)
            assert info.length == len(large_seq)
            
            rev_comp = mcp_server.reverse_complement(large_seq)
            assert len(rev_comp["reverse_complement"]) == len(large_seq)
    
    def test_multiple_operations_memory(self, mcp_server):
        """Test memory usage across multiple operations."""
        sequence = BACTERIAL_GENES["lacZ"]["sequence"]
        
        # Perform multiple operations in sequence
        operations_count = 0
        
        for _ in range(10):  # Repeat operations
            info = mcp_server.sequence_info(sequence)
            rev_comp = mcp_server.reverse_complement(sequence)
            restriction = mcp_server.restriction_analysis(sequence, ["EcoRI", "BamHI"])
            operations_count += 3
        
        # All operations should complete
        assert operations_count == 30
    
    def test_concurrent_operations_simulation(self, mcp_server):
        """Simulate concurrent operations (sequential but rapid)."""
        sequences = [
            BACTERIAL_GENES["lacZ"]["sequence"][:200],
            BACTERIAL_GENES["ampR"]["sequence"][:200],
            PERFORMANCE_TEST_SEQUENCES["gfp_like"][:200]
        ]
        
        start_time = time.time()
        
        # Simulate rapid sequential operations
        results = []
        for seq in sequences:
            for _ in range(3):  # 3 operations per sequence
                r1 = mcp_server.sequence_info(seq)
                r2 = mcp_server.reverse_complement(seq)
                r3 = mcp_server.restriction_analysis(seq, ["EcoRI"])
                results.extend([r1, r2, r3])
        
        end_time = time.time()
        
        # All operations completed
        assert len(results) == len(sequences) * 3 * 3
        # Should complete reasonably quickly
        assert end_time - start_time < 10.0


class TestWorstCaseScenarios:
    """Test worst-case performance scenarios."""
    
    def test_high_gc_content_sequence(self, mcp_server):
        """Test performance with high GC content sequences."""
        high_gc_seq = "GCGCGCGC" * 1000  # 8kb of alternating GC
        
        start_time = time.time()
        info = mcp_server.sequence_info(high_gc_seq)
        end_time = time.time()
        
        assert end_time - start_time < 3.0
        assert info.gc_content > 95.0  # Should be very high GC
    
    def test_low_complexity_sequence(self, mcp_server):
        """Test performance with low complexity sequences."""
        low_complexity = "AAAA" * 2000  # 8kb of just A's
        
        start_time = time.time()
        info = mcp_server.sequence_info(low_complexity)
        rev_comp = mcp_server.reverse_complement(low_complexity)
        end_time = time.time()
        
        assert end_time - start_time < 2.0
        assert info.gc_content < 5.0  # Should be very low GC
        assert rev_comp["reverse_complement"] == "TTTT" * 2000
    
    def test_many_restriction_sites(self, mcp_server):
        """Test performance with sequences containing many restriction sites."""
        # Create sequence with many EcoRI sites
        ecori_rich = ("ATGC" * 5 + "GAATTC") * 500  # ~500 EcoRI sites in ~8kb
        
        start_time = time.time()
        results = mcp_server.restriction_analysis(ecori_rich, ["EcoRI"])
        digest_result = mcp_server.digest_sequence(ecori_rich, ["EcoRI"])
        end_time = time.time()
        
        assert end_time - start_time < 5.0
        # Should find many cutting sites
        ecori_result = [r for r in results if r.enzyme_name == "EcoRI"][0]
        assert len(ecori_result.cut_positions) > 400  # Should find most sites
        assert digest_result["num_fragments"] > 400
    
    def test_primers_with_many_binding_sites(self, mcp_server):
        """Test primer annealing with sequences having many potential binding sites."""
        # Create template with repetitive sequences
        repetitive_template = ("ATGCGATCGCTA" * 100 +  # 1.2kb of repetitive sequence
                              "GAATTCATGCGATCGCTAGAATTC" +  # Unique middle section
                              "ATGCGATCGCTA" * 100)  # Another 1.2kb
        
        primers = ["ATGCGATCGCTA", "TAGATCGCTATG"]  # Primers that match repetitive parts
        
        start_time = time.time()
        annealing = mcp_server.anneal_primers(primers, repetitive_template, limit=10)
        end_time = time.time()
        
        assert end_time - start_time < 10.0
        # Should find multiple binding sites
        assert annealing["products"] >= 1


@pytest.mark.slow
class TestLongRunningOperations:
    """Test operations that are expected to take longer."""
    
    def test_very_large_sequence_analysis(self, mcp_server):
        """Test analysis of very large sequences (like small genomes)."""
        # Simulate a small bacterial genome fragment (~100kb)
        large_genome_fragment = (BACTERIAL_GENES["lacZ"]["sequence"] * 50 +
                               BACTERIAL_GENES["ampR"]["sequence"] * 30)[:100000]
        
        # This is allowed to take longer
        start_time = time.time()
        info = mcp_server.sequence_info(large_genome_fragment)
        end_time = time.time()
        
        # Allow up to 30 seconds for very large sequences
        assert end_time - start_time < 30.0
        assert info.length == len(large_genome_fragment)
        assert 0 <= info.gc_content <= 100
    
    def test_comprehensive_restriction_mapping(self, mcp_server):
        """Test comprehensive restriction mapping with many enzymes."""
        genome_fragment = PERFORMANCE_TEST_SEQUENCES["large_genome_fragment"]
        
        # Use many common restriction enzymes
        many_enzymes = [
            "EcoRI", "BamHI", "HindIII", "XbaI", "SacI", "KpnI", "SmaI", "PstI",
            "NcoI", "NotI", "XhoI", "SpeI", "ApaI", "SalI", "ClaI", "BglII",
            "AvrII", "NheI", "PmeI", "ScaI", "EcoRV", "HpaI", "DraI", "SspI"
        ]
        
        start_time = time.time()
        results = mcp_server.restriction_analysis(genome_fragment, many_enzymes)
        end_time = time.time()
        
        # Allow up to 60 seconds for comprehensive analysis
        assert end_time - start_time < 60.0
        assert len(results) <= len(many_enzymes)
    
    def test_complex_multi_fragment_assembly(self, mcp_server):
        """Test complex assembly with many fragments."""
        # Create 15 fragments for complex assembly
        base_sequences = [
            BACTERIAL_GENES["lacZ"]["sequence"][:200],
            BACTERIAL_GENES["ampR"]["sequence"][:200],
            PERFORMANCE_TEST_SEQUENCES["gfp_like"][:200]
        ]
        
        fragments = []
        overlap = "GGGGCCCCAAAATTTT"  # 16bp overlap
        
        for i in range(15):
            base_seq = base_sequences[i % len(base_sequences)]
            if i == 0:
                fragment = base_seq + overlap
            elif i == 14:
                fragment = overlap + base_seq
            else:
                fragment = overlap + base_seq + overlap
            fragments.append(fragment)
        
        start_time = time.time()
        try:
            result = mcp_server.assembly(fragments, limit=15)
            end_time = time.time()
            
            # Allow up to 120 seconds for complex assembly
            assert end_time - start_time < 120.0
            assert isinstance(result, AssemblyResult)
        except ValueError:
            # Complex assembly might fail - that's acceptable
            end_time = time.time()
            # Should fail reasonably quickly, not hang
            assert end_time - start_time < 30.0


@pytest.mark.benchmark
class TestBenchmarkComparisons:
    """Benchmark tests for comparing performance across versions."""
    
    def test_sequence_operations_benchmark(self, mcp_server):
        """Benchmark basic sequence operations."""
        test_sequence = BACTERIAL_GENES["lacZ"]["sequence"]
        
        operations = []
        
        # Time each operation
        for operation_name, operation in [
            ("sequence_info", lambda: mcp_server.sequence_info(test_sequence)),
            ("reverse_complement", lambda: mcp_server.reverse_complement(test_sequence)),
            ("restriction_analysis", lambda: mcp_server.restriction_analysis(test_sequence, ["EcoRI", "BamHI"])),
        ]:
            start_time = time.time()
            result = operation()
            end_time = time.time()
            
            operations.append({
                "operation": operation_name,
                "time": end_time - start_time,
                "result_type": type(result).__name__
            })
        
        # Verify all operations completed
        assert len(operations) == 3
        
        # Print benchmark results (will appear in pytest output with -s flag)
        print("\nBenchmark Results:")
        for op in operations:
            print(f"  {op['operation']}: {op['time']:.4f}s -> {op['result_type']}")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-m", "not slow"]) 