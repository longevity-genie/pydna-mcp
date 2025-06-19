#!/usr/bin/env python3
"""Test script to verify real pydna functionality for the MCP server."""

# Test the basic pydna functionality that we'll expose through MCP
import sys
import os

from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydna.design import primer_design
from pydna.amplify import pcr, Anneal
from pydna.assembly import Assembly
from pydna.tm import tm_default
from Bio.Restriction import EcoRI, BamHI, XhoI

def test_basic_sequence():
    """Test basic sequence creation and properties."""
    print("=== Testing Basic Sequence ===")
    
    # Create a test sequence
    seq_str = "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT"
    dseq = Dseqrecord(seq_str, name="test_seq")
    
    print(f"Sequence: {dseq.seq}")
    print(f"Length: {len(dseq)}")
    print(f"Circular: {dseq.circular}")
    print(f"Molecular weight: {dseq.seq.mw()}")
    
    # GC content calculation
    gc_content = (seq_str.count('G') + seq_str.count('C')) / len(seq_str) * 100
    print(f"GC content: {gc_content:.2f}%")
    
    # Melting temperature
    tm = tm_default(seq_str)
    print(f"Tm: {tm:.2f}°C")
    
    # Reverse complement
    rc = dseq.reverse_complement()
    print(f"Reverse complement: {rc.seq}")
    
    return dseq

def test_primer_design():
    """Test primer design functionality."""
    print("\n=== Testing Primer Design ===")
    
    # Use a longer template that's more suitable for primer design
    template = "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAATAGGATCCTTTTTTTTAAAA"
    template_seq = Dseqrecord(template)
    
    try:
        # Design primers
        amplicon = primer_design(template_seq, target_tm=55, limit=13)
        
        print(f"Template: {template}")
        print(f"Forward primer: {amplicon.forward_primer.seq}")
        print(f"Reverse primer: {amplicon.reverse_primer.seq}")
        print(f"Amplicon length: {len(amplicon)}")
        
        # Calculate primer properties
        fp = str(amplicon.forward_primer.seq)
        rp = str(amplicon.reverse_primer.seq)
        
        fp_tm = tm_default(fp)
        rp_tm = tm_default(rp)
        
        print(f"Forward primer Tm: {fp_tm:.2f}°C")
        print(f"Reverse primer Tm: {rp_tm:.2f}°C")
        
        return amplicon
    except Exception as e:
        print(f"Primer design failed: {e}")
        return None

def test_pcr():
    """Test PCR simulation."""
    print("\n=== Testing PCR ===")
    
    template = "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAATAGGATCCTTTTTTTTAAAA"
    template_seq = Dseqrecord(template)
    
    try:
        # Use manually designed primers that we know will work
        fp = "ATGAAAGCACTGATTC"  # Forward primer
        rp = "TTTTAAAAAAAAGGA"    # Reverse primer (reverse complement of end)
        
        # Perform PCR
        amplicon = pcr(fp, rp, template_seq)
        
        print(f"Template: {template}")
        print(f"Forward primer: {fp}")
        print(f"Reverse primer: {rp}")
        print(f"PCR product: {amplicon.seq}")
        print(f"Product length: {len(amplicon)}")
        
        return amplicon
    except Exception as e:
        print(f"PCR failed: {e}")
        return None

def test_restriction():
    """Test restriction enzyme analysis."""
    print("\n=== Testing Restriction Analysis ===")
    
    # Use a sequence with known restriction sites
    sequence = "GAATTCATGAAAGCACTGATTCTATTGCTGAAAAAGATAATAGGATCC"
    dseq = Dseqrecord(sequence)
    
    print(f"Sequence: {sequence}")
    
    # Test with EcoRI (GAATTC)
    ecori_sites = EcoRI.search(dseq.seq)
    print(f"EcoRI sites: {ecori_sites}")
    
    # Test with BamHI (GGATCC)  
    bamhi_sites = BamHI.search(dseq.seq)
    print(f"BamHI sites: {bamhi_sites}")
    
    # Digest with EcoRI
    if ecori_sites:
        fragments = dseq.cut(EcoRI)
        print(f"EcoRI digestion fragments: {[len(f) for f in fragments]}")
        for i, frag in enumerate(fragments):
            print(f"  Fragment {i+1}: {frag.seq} (length: {len(frag)})")
    
    # Digest with BamHI
    if bamhi_sites:
        fragments = dseq.cut(BamHI)
        print(f"BamHI digestion fragments: {[len(f) for f in fragments]}")
        for i, frag in enumerate(fragments):
            print(f"  Fragment {i+1}: {frag.seq} (length: {len(frag)})")
    
    return dseq

def test_assembly():
    """Test DNA assembly."""
    print("\n=== Testing Assembly ===")
    
    # Create overlapping fragments for assembly
    frag1 = "ATGAAAGCACTGATTCTATTGC"        # 22 bp
    frag2 = "CTATTGCTGAAAAAGATAAT"          # 20 bp, overlaps 8 bp with frag1
    frag3 = "AGATAATAGGATCCTTTTTTTT"        # 22 bp, overlaps 8 bp with frag2
    
    fragments = [Dseqrecord(frag1), Dseqrecord(frag2), Dseqrecord(frag3)]
    
    print(f"Fragment 1: {frag1}")
    print(f"Fragment 2: {frag2}")  
    print(f"Fragment 3: {frag3}")
    
    try:
        # Perform assembly with smaller overlap limit
        assembly = Assembly(fragments, limit=8)
        linear_products = assembly.assemble_linear()
        
        if linear_products:
            product = linear_products[0]
            print(f"Assembly product: {product.seq}")
            print(f"Product length: {len(product)}")
            print(f"Is circular: {product.circular}")
        else:
            print("No linear assembly products found")
            # Try circular assembly
            circular_products = assembly.assemble_circular()
            if circular_products:
                product = circular_products[0]
                print(f"Circular assembly product: {product.seq}")
                print(f"Product length: {len(product)}")
                print(f"Is circular: {product.circular}")
            else:
                print("No assembly products found")
        
        return assembly
    except Exception as e:
        print(f"Assembly failed: {e}")
        return None

def test_annealing():
    """Test primer annealing analysis."""
    print("\n=== Testing Primer Annealing ===")
    
    template = "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT"
    template_seq = Dseqrecord(template)
    
    # Test primers - use primers that should anneal
    primers = [
        Dseqrecord("ATGAAAGCACTGATT"),         # Forward primer (first 15 bp)
        Dseqrecord("ATTATCTTTTTCAGC"),         # This is the reverse complement of the end
    ]
    
    try:
        annealing = Anneal(primers, template_seq, limit=10)
        
        print(f"Template: {template}")
        print(f"Forward primers found: {len(annealing.forward_primers)}")
        print(f"Reverse primers found: {len(annealing.reverse_primers)}")
        
        for fp in annealing.forward_primers:
            print(f"  Forward primer at position {fp.position}: {fp.seq}")
        
        for rp in annealing.reverse_primers:
            print(f"  Reverse primer at position {rp.position}: {rp.seq}")
        
        return annealing
    except Exception as e:
        print(f"Annealing analysis failed: {e}")
        return None

def test_simple_pcr():
    """Test simple PCR with basic sequences."""
    print("\n=== Testing Simple PCR ===")
    
    # Simple template and primers for basic PCR test
    template = "ATGCGATCGTAGCTAGCTAGCTAGCATGC"
    template_seq = Dseqrecord(template)
    
    # Simple primers
    forward = "ATGCGATCGTAGC"
    reverse = "GCATGCTAGCTAGC"  # This should be the reverse complement 
    
    try:
        # Try PCR
        result = pcr(forward, reverse, template_seq)
        print(f"PCR successful!")
        print(f"Product: {result.seq}")
        print(f"Length: {len(result)}")
        return result
    except Exception as e:
        print(f"Simple PCR failed: {e}")
        return None

def main():
    """Run all PyDNA functionality tests."""
    print("Starting PyDNA functionality tests...")
    
    # Run tests
    seq_result = test_basic_sequence()
    primer_result = test_primer_design()
    pcr_result = test_pcr()
    restriction_result = test_restriction()
    assembly_result = test_assembly()
    annealing_result = test_annealing()
    simple_pcr_result = test_simple_pcr()
    
    print("\n=== Test Summary ===")
    print(f"Basic sequence: {'✓' if seq_result else '✗'}")
    print(f"Primer design: {'✓' if primer_result else '✗'}")
    print(f"PCR simulation: {'✓' if pcr_result else '✗'}")
    print(f"Restriction analysis: {'✓' if restriction_result else '✗'}")
    print(f"Assembly: {'✓' if assembly_result else '✗'}")
    print(f"Annealing: {'✓' if annealing_result else '✗'}")
    print(f"Simple PCR: {'✓' if simple_pcr_result else '✗'}")

if __name__ == "__main__":
    main() 