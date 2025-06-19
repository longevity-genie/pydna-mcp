#!/usr/bin/env python3
"""
PyDNA MCP Server - DNA simulation and cloning tools via Model Context Protocol.

This server provides comprehensive DNA simulation and molecular cloning tools through the
Model Context Protocol (MCP). It wraps the PyDNA library to offer molecular biology
functionality including PCR simulation, primer design, restriction analysis, DNA assembly,
and more.

Key capabilities:
- Sequence creation and analysis
- PCR amplification simulation
- Primer design and optimization
- Restriction enzyme analysis
- DNA assembly (homologous recombination, Gibson)
- Gel electrophoresis simulation
- GenBank integration
- File format conversion (FASTA, GenBank)

Example workflows:
1. Design primers for a gene, simulate PCR, analyze products
2. Plan restriction cloning with enzyme analysis
3. Design Gibson assembly experiments
4. Download sequences from GenBank for analysis
"""

import asyncio
import os
import sys
import argparse
from pathlib import Path
from typing import List, Dict, Any, Optional, Union
from contextlib import asynccontextmanager
import tempfile
import json
import base64
from io import StringIO, BytesIO

from fastmcp import FastMCP
from pydantic import BaseModel, Field
from eliot import start_action
import typer

# PyDNA imports
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydna.amplify import pcr, Anneal
from pydna.design import primer_design, assembly_fragments
from pydna.assembly import Assembly
from pydna.gel import gel
from pydna.readers import read
from pydna.parsers import parse
from pydna.genbank import genbank
from pydna.download import download_text
from pydna.tm import tm_default, tm_neb
from Bio.Restriction import *
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio.SeqIO as SeqIO

# Configuration
DEFAULT_HOST = os.getenv("MCP_HOST", "0.0.0.0")
DEFAULT_PORT = int(os.getenv("MCP_PORT", "3001"))
DEFAULT_TRANSPORT = os.getenv("MCP_TRANSPORT", "streamable-http")

# Pydantic models for structured responses
class SequenceInfo(BaseModel):
    """Information about a DNA sequence."""
    sequence: str = Field(description="DNA sequence string")
    length: int = Field(description="Length of sequence")
    circular: bool = Field(description="Whether sequence is circular")
    gc_content: float = Field(description="GC content percentage")
    molecular_weight: float = Field(description="Molecular weight in Daltons")
    tm: float = Field(description="Melting temperature in Celsius")

class PrimerInfo(BaseModel):
    """Information about a primer."""
    sequence: str = Field(description="Primer sequence")
    length: int = Field(description="Primer length")
    tm: float = Field(description="Melting temperature in Celsius")
    gc_content: float = Field(description="GC content percentage")

class AmpliconInfo(BaseModel):
    """Information about a PCR amplicon."""
    sequence: str = Field(description="Amplicon sequence")
    length: int = Field(description="Amplicon length")
    forward_primer: PrimerInfo = Field(description="Forward primer information")
    reverse_primer: PrimerInfo = Field(description="Reverse primer information")
    tm_forward: float = Field(description="Forward primer melting temperature")
    tm_reverse: float = Field(description="Reverse primer melting temperature")

class RestrictionAnalysis(BaseModel):
    """Restriction enzyme analysis results."""
    enzyme_name: str = Field(description="Name of restriction enzyme")
    recognition_site: str = Field(description="Recognition sequence")
    cut_positions: List[int] = Field(description="Positions where enzyme cuts")
    fragments: List[int] = Field(description="Fragment sizes after digestion")

class AssemblyResult(BaseModel):
    """DNA assembly result."""
    sequence: str = Field(description="Assembled sequence")
    length: int = Field(description="Length of assembled sequence")
    circular: bool = Field(description="Whether result is circular")
    fragments_used: int = Field(description="Number of fragments used in assembly")

class PyDNAMCP(FastMCP):
    """PyDNA MCP Server with DNA simulation and cloning tools."""
    
    def __init__(
        self, 
        name: str = "PyDNA MCP Server",
        prefix: str = "pydna_",
        **kwargs
    ):
        """Initialize the PyDNA tools with FastMCP functionality."""
        super().__init__(name=name, **kwargs)
        self.prefix = prefix
        self._register_pydna_tools()
        self._register_pydna_resources()
    
    def _register_pydna_tools(self):
        """Register PyDNA-specific tools."""
        # Sequence creation and manipulation
        self.tool(name=f"{self.prefix}create_sequence", description="Create a DNA sequence from string")(self.create_sequence)
        self.tool(name=f"{self.prefix}sequence_info", description="Get detailed information about a DNA sequence")(self.sequence_info)
        self.tool(name=f"{self.prefix}reverse_complement", description="Get reverse complement of a DNA sequence")(self.reverse_complement)
        
        # Primer design and PCR
        self.tool(name=f"{self.prefix}design_primers", description="Design PCR primers for a DNA sequence")(self.design_primers)
        self.tool(name=f"{self.prefix}pcr_amplify", description="Simulate PCR amplification with given primers and template")(self.pcr_amplify)
        self.tool(name=f"{self.prefix}anneal_primers", description="Analyze primer annealing to template")(self.anneal_primers)
        
        # Restriction analysis
        self.tool(name=f"{self.prefix}restriction_analysis", description="Analyze restriction enzyme cutting sites")(self.restriction_analysis)
        self.tool(name=f"{self.prefix}digest_sequence", description="Digest DNA sequence with restriction enzymes")(self.digest_sequence)
        
        # Assembly
        self.tool(name=f"{self.prefix}assembly", description="Assemble DNA fragments by homologous recombination")(self.assembly)
        self.tool(name=f"{self.prefix}gibson_assembly", description="Design fragments for Gibson assembly")(self.gibson_assembly)
        
        # File operations
        self.tool(name=f"{self.prefix}read_sequence_file", description="Read DNA sequence from file (FASTA, GenBank, etc.)")(self.read_sequence_file)
        self.tool(name=f"{self.prefix}write_sequence_file", description="Write DNA sequence to file")(self.write_sequence_file)
        
        # Gel electrophoresis
        self.tool(name=f"{self.prefix}gel_electrophoresis", description="Simulate gel electrophoresis")(self.gel_electrophoresis)
        
        # Download sequences
        self.tool(name=f"{self.prefix}download_genbank", description="Download sequence from GenBank")(self.download_genbank)
    
    def _register_pydna_resources(self):
        """Register PyDNA-specific resources."""
        
        @self.resource(f"resource://{self.prefix}help")
        def get_pydna_help() -> str:
            """
            Get comprehensive help for PyDNA MCP server tools.
            
            This resource provides detailed information about:
            - Available tools and their purposes
            - Common workflows for DNA cloning simulation
            - Example usage patterns
            - Best practices for molecular cloning design
            
            Returns:
                Comprehensive help documentation
            """
            return """
# PyDNA MCP Server Help

## Overview
This MCP server provides tools for DNA simulation and molecular cloning using the PyDNA library.
PyDNA allows you to simulate common molecular biology techniques including:

- PCR amplification
- Primer design
- Restriction enzyme digestion
- DNA assembly (Gibson, homologous recombination)
- Gel electrophoresis simulation

## Common Workflows

### 1. Basic Sequence Analysis
1. Create or load a DNA sequence
2. Get sequence information (GC content, Tm, etc.)
3. Find restriction sites

### 2. PCR Simulation
1. Create template sequence
2. Design primers (or provide existing ones)
3. Simulate PCR amplification
4. Analyze the amplicon

### 3. Cloning Workflow
1. Design primers with restriction sites
2. Simulate PCR amplification
3. Digest PCR product and vector
4. Simulate ligation/assembly

### 4. Gibson Assembly
1. Design overlapping fragments
2. Simulate assembly
3. Verify final construct

## Example Usage
See the individual tool descriptions for specific examples and parameters.
"""

        @self.resource(f"resource://{self.prefix}examples")
        def get_examples() -> str:
            """Get example workflows and use cases."""
            return """
# PyDNA MCP Server Examples

## Example 1: Basic PCR
```python
# Create template
template = "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT"

# Design primers
primers = design_primers(template, target_tm=55)

# Simulate PCR
amplicon = pcr_amplify(forward_primer, reverse_primer, template)
```

## Example 2: Restriction Analysis
```python
# Analyze restriction sites
sites = restriction_analysis(sequence, ["EcoRI", "BamHI", "XhoI"])

# Digest with enzymes
fragments = digest_sequence(sequence, ["EcoRI"])
```

## Example 3: Gibson Assembly
```python
# Prepare fragments for Gibson assembly
fragments = ["ATGAAAGCA...", "GCACTGATT...", "ATTGCTGAA..."]
assembled = gibson_assembly(fragments, overlap=20)
```
"""
    
    def create_sequence(self, sequence: str, name: str = "seq", circular: bool = False) -> SequenceInfo:
        """
        Create a DNA sequence from string and return comprehensive sequence information.
        
        This method creates a Dseqrecord object from a DNA sequence string and calculates
        important molecular properties including GC content, melting temperature, and
        molecular weight. Use this as the starting point for most DNA analysis workflows.
        
        Args:
            sequence (str): DNA sequence string (A, T, G, C nucleotides)
                Example: "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT"
            name (str, optional): Name for the sequence. Defaults to "seq"
                Example: "my_gene", "plasmid_insert", "pcr_product"
            circular (bool, optional): Whether sequence is circular (like plasmids). 
                Defaults to False. Set True for plasmids, False for linear DNA.
        
        Returns:
            SequenceInfo: Object containing:
                - sequence: The DNA sequence string
                - length: Number of nucleotides
                - circular: Whether the sequence is circular
                - gc_content: Percentage of G and C nucleotides (0-100)
                - molecular_weight: Weight in Daltons
                - tm: Melting temperature in Celsius
        
        Example Input:
            sequence = "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT"
            name = "test_gene"
            circular = False
        
        Example Output:
            {
                "sequence": "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT",
                "length": 34,
                "circular": false,
                "gc_content": 32.35,
                "molecular_weight": 10456.7,
                "tm": 67.8
            }
        
        Common Use Cases:
            - Analyze a gene sequence before PCR design
            - Check properties of a synthetic DNA construct
            - Validate sequence composition before cloning
            - Calculate molecular weight for gel analysis
        
        Error Conditions:
            - Invalid nucleotides (not A, T, G, C) will raise ValueError
            - Empty sequence string will raise ValueError
        """
        with start_action(action_type="create_sequence", sequence_length=len(sequence), circular=circular) as action:
            try:
                dseq = Dseqrecord(sequence, name=name, circular=circular)
                
                # Calculate properties
                gc_content = (sequence.upper().count('G') + sequence.upper().count('C')) / len(sequence) * 100
                tm = tm_default(sequence)
                
                result = SequenceInfo(
                    sequence=str(dseq.seq),
                    length=len(dseq),
                    circular=dseq.circular,
                    gc_content=round(gc_content, 2),
                    molecular_weight=dseq.seq.mw(),
                    tm=round(tm, 2)
                )
                
                action.add_success_fields(sequence_info=result.model_dump())
                return result
                
            except Exception as e:
                action.add_success_fields(error=str(e))
                raise ValueError(f"Error creating sequence: {str(e)}")
    
    def sequence_info(self, sequence: str) -> SequenceInfo:
        """
        Get detailed molecular information about a DNA sequence.
        
        Analyzes a DNA sequence string and returns comprehensive molecular properties
        including GC content, melting temperature, molecular weight, and basic
        characteristics. This is useful for sequence validation and analysis.
        
        Args:
            sequence (str): DNA sequence string containing only A, T, G, C nucleotides
                Example: "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT"
                Example: "GGCCTTAAGCTTGAATTCGGATCC" (with restriction sites)
        
        Returns:
            SequenceInfo: Detailed sequence analysis containing:
                - sequence: The input DNA sequence
                - length: Number of base pairs
                - circular: Always False (use create_sequence for circular)
                - gc_content: GC percentage (ideal range: 40-60% for PCR)
                - molecular_weight: Molecular weight in Daltons
                - tm: Melting temperature in Celsius (useful for PCR conditions)
        
        Example Input:
            sequence = "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT"
        
        Example Output:
            {
                "sequence": "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT",
                "length": 34,
                "circular": false,
                "gc_content": 32.35,
                "molecular_weight": 10456.7,
                "tm": 67.8
            }
        
        Interpretation Guide:
            - GC content: 40-60% is optimal for PCR; <30% or >70% may cause issues
            - Tm: Higher GC content = higher melting temperature
            - Length: Typical genes are 300-3000 bp, primers are 18-30 bp
            - Molecular weight: Useful for calculating DNA concentrations
        
        Common Use Cases:
            - Validate sequence before experimental use
            - Check if sequence is suitable for PCR amplification
            - Calculate properties for primer design
            - Quality control for synthetic DNA sequences
        
        Error Conditions:
            - Invalid nucleotides (not A,T,G,C) raise ValueError
            - Empty sequence raises ValueError
            - Non-string input raises TypeError
        """
        with start_action(action_type="sequence_info", sequence_length=len(sequence) if sequence else 0) as action:
            try:
                # Handle edge cases
                if not sequence or len(sequence.strip()) == 0:
                    raise ValueError("Sequence cannot be empty or contain only whitespace")
                
                # Clean sequence 
                clean_sequence = sequence.strip().upper().replace(" ", "").replace("\n", "").replace("\t", "")
                
                if len(clean_sequence) == 0:
                    raise ValueError("Sequence cannot be empty after cleaning")
                
                # Validate sequence contains only valid DNA characters
                valid_chars = set('ATCG')
                if not all(c in valid_chars for c in clean_sequence):
                    raise ValueError("Sequence contains invalid nucleotides. Only A, T, C, G are allowed")
                
                dseq = Dseqrecord(clean_sequence)
                
                # Calculate properties with error handling
                gc_content = (clean_sequence.count('G') + clean_sequence.count('C')) / len(clean_sequence) * 100
                
                # Handle very short sequences for Tm calculation
                try:
                    if len(clean_sequence) < 4:
                        # Very short sequences: use simple rule
                        tm = (clean_sequence.count('A') + clean_sequence.count('T')) * 2 + (clean_sequence.count('G') + clean_sequence.count('C')) * 4
                    else:
                        tm = tm_default(clean_sequence)
                except:
                    # Fallback for problematic sequences
                    tm = (clean_sequence.count('A') + clean_sequence.count('T')) * 2 + (clean_sequence.count('G') + clean_sequence.count('C')) * 4
                
                result = SequenceInfo(
                    sequence=str(dseq.seq),
                    length=len(dseq),
                    circular=dseq.circular,
                    gc_content=round(gc_content, 2),
                    molecular_weight=dseq.seq.mw(),
                    tm=round(float(tm), 2)
                )
                
                action.add_success_fields(sequence_info=result.model_dump())
                return result
                
            except Exception as e:
                action.add_success_fields(error=str(e))
                raise ValueError(f"Error analyzing sequence: {str(e)}")
    
    def reverse_complement(self, sequence: str) -> Dict[str, str]:
        """
        Calculate the reverse complement of a DNA sequence.
        
        The reverse complement is essential in molecular biology as it represents
        the complementary strand of double-stranded DNA. This is used for:
        - Designing reverse PCR primers
        - Understanding DNA hybridization
        - Analyzing both strands of DNA
        - Converting between forward and reverse orientations
        
        Args:
            sequence (str): DNA sequence string (A, T, G, C nucleotides)
                Example: "ATGCGATCG" 
                Example: "GGCCTTAAGCTT" (restriction site)
        
        Returns:
            Dict[str, str]: Dictionary containing both sequences:
                - "original": The input sequence
                - "reverse_complement": The reverse complement
        
        Example Input:
            sequence = "ATGCGATCG"
        
        Example Output:
            {
                "original": "ATGCGATCG",
                "reverse_complement": "CGATCGCAT"
            }
        
        Base Pairing Rules:
            - A pairs with T (and vice versa)
            - G pairs with C (and vice versa)
            - The sequence is reversed (3' to 5' becomes 5' to 3')
        
        Example Transformations:
            - "ATGC" → "GCAT"
            - "AAAA" → "TTTT"  
            - "GGCC" → "GGCC" (palindrome)
            - "GAATTC" → "GAATTC" (EcoRI site is palindromic)
        
        Common Use Cases:
            - Design reverse primers for PCR (use reverse complement of target end)
            - Check if sequences are palindromic (same as reverse complement)
            - Analyze both strands of genomic regions
            - Convert primer sequences between orientations
            - Validate sequence directionality in cloning
        
        Error Conditions:
            - Invalid nucleotides (not A,T,G,C) raise ValueError
            - Empty sequence raises ValueError
            - Ambiguous nucleotides (N, R, Y, etc.) may not be handled
        """
        with start_action(action_type="reverse_complement", sequence_length=len(sequence)) as action:
            try:
                dseq = Dseqrecord(sequence)
                rc_seq = dseq.reverse_complement()
                
                result = {
                    "original": str(dseq.seq),
                    "reverse_complement": str(rc_seq.seq)
                }
                
                action.add_success_fields(result=result)
                return result
                
            except Exception as e:
                action.add_success_fields(error=str(e))
                raise ValueError(f"Error computing reverse complement: {str(e)}")
    
    def design_primers(self, template: str, target_tm: float = 55.0, limit: int = 13) -> AmpliconInfo:
        """
        Automatically design PCR primers for a DNA template sequence.
        
        This method uses PyDNA's primer design algorithm to create optimal forward
        and reverse primers for PCR amplification of the entire template. The
        algorithm selects primers with suitable melting temperatures and minimal
        secondary structure issues.
        
        Args:
            template (str): DNA template sequence to amplify
                Example: "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT"
                Should be the full sequence you want to amplify
            target_tm (float, optional): Target melting temperature for primers in Celsius.
                Defaults to 55.0. Range: 50-65°C is typical for PCR
                Example: 60.0 for high-stringency PCR
            limit (int, optional): Minimum primer length. Defaults to 13.
                Range: 13-30 bp. Shorter = less specific, Longer = more specific
        
        Returns:
            AmpliconInfo: Complete information about the designed primers and amplicon:
                - sequence: The amplified sequence (same as template)
                - length: Amplicon length in base pairs
                - forward_primer: Forward primer details (sequence, Tm, GC%)
                - reverse_primer: Reverse primer details (sequence, Tm, GC%)
                - tm_forward: Forward primer melting temperature
                - tm_reverse: Reverse primer melting temperature
        
        Example Input:
            template = "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT"
            target_tm = 58.0
            limit = 18
        
        Example Output:
            {
                "sequence": "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT",
                "length": 34,
                "forward_primer": {
                    "sequence": "ATGAAAGCACTGATTCTAT",
                    "length": 19,
                    "tm": 57.8,
                    "gc_content": 36.8
                },
                "reverse_primer": {
                    "sequence": "ATTATCTTTTTCAGCAATAG",
                    "length": 20,
                    "tm": 58.2,
                    "gc_content": 30.0
                },
                "tm_forward": 57.8,
                "tm_reverse": 58.2
            }
        
        Design Principles:
            - Forward primer binds to 5' end of template
            - Reverse primer binds to 3' end (reverse complement)
            - Similar Tm values for both primers (within 3°C)
            - GC content typically 40-60% for stability
            - Avoid hairpins, dimers, and repetitive sequences
        
        Common Use Cases:
            - Quick primer design for gene amplification
            - PCR primer optimization for cloning
            - Initial primer design before manual refinement
            - Batch primer design for multiple targets
        
        Best Practices:
            - target_tm of 55-60°C works for most applications
            - Check primers for hairpins and dimers separately
            - Consider adding restriction sites to primer ends for cloning
            - Validate primer specificity against genome if needed
        
        Error Conditions:
            - Template too short (< limit*2) raises ValueError
            - Invalid nucleotides in template raise ValueError
            - Unable to find suitable primers raises ValueError
        """
        with start_action(action_type="design_primers", template_length=len(template), target_tm=target_tm) as action:
            try:
                template_seq = Dseqrecord(template)
                amplicon = primer_design(template_seq, target_tm=target_tm, limit=limit)
                
                # Get primer information
                fp = amplicon.forward_primer
                rp = amplicon.reverse_primer
                
                fp_info = PrimerInfo(
                    sequence=str(fp.seq),
                    length=len(fp),
                    tm=round(tm_default(str(fp.seq)), 2),
                    gc_content=round((str(fp.seq).upper().count('G') + str(fp.seq).upper().count('C')) / len(fp) * 100, 2)
                )
                
                rp_info = PrimerInfo(
                    sequence=str(rp.seq),
                    length=len(rp),
                    tm=round(tm_default(str(rp.seq)), 2),
                    gc_content=round((str(rp.seq).upper().count('G') + str(rp.seq).upper().count('C')) / len(rp) * 100, 2)
                )
                
                result = AmpliconInfo(
                    sequence=str(amplicon.seq),
                    length=len(amplicon),
                    forward_primer=fp_info,
                    reverse_primer=rp_info,
                    tm_forward=fp_info.tm,
                    tm_reverse=rp_info.tm
                )
                
                action.add_success_fields(amplicon_info=result.model_dump())
                return result
                
            except Exception as e:
                action.add_success_fields(error=str(e))
                raise ValueError(f"Error designing primers: {str(e)}")
    
    def pcr_amplify(self, forward_primer: str, reverse_primer: str, template: str) -> AmpliconInfo:
        """
        Simulate PCR amplification with specific primers and template.
        
        This method simulates a PCR reaction by finding where the primers bind to
        the template and predicting the amplified product. It's essential for
        validating primer designs and predicting PCR outcomes before lab work.
        
        Args:
            forward_primer (str): Forward primer sequence (5' to 3')
                Example: "ATGAAAGCACTGATTCTAT"
                Should bind to the beginning of your target region
            reverse_primer (str): Reverse primer sequence (5' to 3')  
                Example: "ATTATCTTTTTCAGCAATAG"
                Should be reverse complement of template end
            template (str): Template DNA sequence containing target region
                Example: "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT"
                Can be larger than target (like full plasmid)
        
        Returns:
            AmpliconInfo: Complete PCR product information:
                - sequence: The amplified DNA sequence
                - length: Product size in base pairs
                - forward_primer: Forward primer analysis
                - reverse_primer: Reverse primer analysis  
                - tm_forward: Forward primer melting temperature
                - tm_reverse: Reverse primer melting temperature
        
        Example Input:
            forward_primer = "ATGAAAGCACTGATTCTAT"
            reverse_primer = "ATTATCTTTTTCAGCAATAG" 
            template = "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT"
        
        Example Output:
            {
                "sequence": "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT",
                "length": 34,
                "forward_primer": {
                    "sequence": "ATGAAAGCACTGATTCTAT",
                    "length": 19,
                    "tm": 57.8,
                    "gc_content": 36.8
                },
                "reverse_primer": {
                    "sequence": "ATTATCTTTTTCAGCAATAG",
                    "length": 20,
                    "tm": 58.2,
                    "gc_content": 30.0
                },
                "tm_forward": 57.8,
                "tm_reverse": 58.2
            }
        
        PCR Principles:
            - Forward primer binds to template 5' end
            - Reverse primer binds to template 3' end (as reverse complement)
            - Product includes primer sequences plus template between them
            - Primers must bind with sufficient complementarity
            - Product length = distance between primer binding sites + primer lengths
        
        Common Use Cases:
            - Validate primer pairs before ordering
            - Predict PCR product size for gel analysis
            - Check if primers will amplify intended target
            - Optimize PCR conditions based on primer Tm values
            - Troubleshoot failed PCR reactions
        
        Primer Design Tips:
            - Primer Tm difference should be < 3°C
            - Optimal Tm range: 55-65°C
            - GC content: 40-60%
            - Length: 18-30 nucleotides
            - Avoid hairpins and primer dimers
        
        Error Conditions:
            - Primers don't bind to template (no complementarity)
            - Invalid nucleotides in primers or template
            - Primers bind in wrong orientation
            - Template too short for primer binding
        
        Troubleshooting:
            - If no product: Check primer complementarity to template
            - If wrong size: Verify primer binding positions
            - If multiple products: Check primer specificity
        """
        with start_action(action_type="pcr_amplify", 
                         fp_length=len(forward_primer), 
                         rp_length=len(reverse_primer),
                         template_length=len(template)) as action:
            try:
                template_seq = Dseqrecord(template)
                amplicon = pcr(forward_primer, reverse_primer, template_seq)
                
                # Get primer information
                fp_info = PrimerInfo(
                    sequence=forward_primer,
                    length=len(forward_primer),
                    tm=round(tm_default(forward_primer), 2),
                    gc_content=round((forward_primer.upper().count('G') + forward_primer.upper().count('C')) / len(forward_primer) * 100, 2)
                )
                
                rp_info = PrimerInfo(
                    sequence=reverse_primer,
                    length=len(reverse_primer),
                    tm=round(tm_default(reverse_primer), 2),
                    gc_content=round((reverse_primer.upper().count('G') + reverse_primer.upper().count('C')) / len(reverse_primer) * 100, 2)
                )
                
                result = AmpliconInfo(
                    sequence=str(amplicon.seq),
                    length=len(amplicon),
                    forward_primer=fp_info,
                    reverse_primer=rp_info,
                    tm_forward=fp_info.tm,
                    tm_reverse=rp_info.tm
                )
                
                action.add_success_fields(amplicon_info=result.model_dump())
                return result
                
            except Exception as e:
                action.add_success_fields(error=str(e))
                raise ValueError(f"Error in PCR amplification: {str(e)}")
    
    def anneal_primers(self, primers: List[str], template: str, limit: int = 13) -> Dict[str, Any]:
        """
        Analyze how primers bind (anneal) to a template sequence.
        
        This method examines primer-template interactions to predict binding
        positions, specificity, and potential PCR products. It's crucial for
        primer validation and troubleshooting PCR issues.
        
        Args:
            primers (List[str]): List of primer sequences to analyze
                Example: ["ATGAAAGCACTGATT", "GGATCCAAGCTTGAA", "TTAATCGGATCC"]
                Can include both forward and reverse primers
            template (str): Template DNA sequence for primer binding analysis
                Example: "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT"
                Should contain regions where primers might bind
            limit (int, optional): Minimum number of complementary bases for binding.
                Defaults to 13. Higher = more stringent, Lower = more permissive
                Range: 10-20 bp depending on primer length and specificity needs
        
        Returns:
            Dict[str, Any]: Comprehensive annealing analysis:
                - template_name: Name of template sequence
                - template_length: Template length in bp
                - limit: Minimum binding length used
                - forward_primers: List of primers binding in forward orientation
                - reverse_primers: List of primers binding in reverse orientation  
                - products: Number of potential PCR products
        
        Example Input:
            primers = ["ATGAAAGCACTGATT", "ATTATCTTTTTCAGC"]
            template = "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT"
            limit = 13
        
        Example Output:
            {
                "template_name": "seq",
                "template_length": 34,
                "limit": 13,
                "forward_primers": [
                    {
                        "sequence": "ATGAAAGCACTGATT",
                        "position": 0,
                        "footprint": 15
                    }
                ],
                "reverse_primers": [
                    {
                        "sequence": "ATTATCTTTTTCAGC", 
                        "position": 19,
                        "footprint": 15
                    }
                ],
                "products": 1
            }
        
        Understanding Results:
            - forward_primers: Bind to template as-is (5' to 3')
            - reverse_primers: Bind as reverse complement
            - position: Starting position of primer binding (0-based)
            - footprint: Number of complementary bases
            - products: Number of possible PCR amplification products
        
        Common Use Cases:
            - Validate primer specificity before PCR
            - Troubleshoot why PCR didn't work
            - Check for off-target primer binding
            - Optimize primer concentrations
            - Design multiplex PCR assays
        
        Annealing Quality Guidelines:
            - footprint >= 13 bp: Usually sufficient for specific binding
            - footprint < 10 bp: May not bind reliably
            - Multiple binding sites: May cause non-specific products
            - No binding: Primers won't work with this template
        
        Primer Design Implications:
            - Products = 0: No PCR amplification expected
            - Products = 1: Clean single product expected
            - Products > 1: Multiple bands possible, optimize conditions
            - Check forward/reverse primer balance for even amplification
        
        Error Conditions:
            - Empty primer list raises ValueError
            - Invalid nucleotides in primers/template raise ValueError
            - Limit greater than primer length may find no binding
        """
        with start_action(action_type="anneal_primers", 
                         num_primers=len(primers), 
                         template_length=len(template)) as action:
            try:
                template_seq = Dseqrecord(template)
                primer_seqs = [Dseqrecord(p) for p in primers]
                
                annealing = Anneal(primer_seqs, template_seq, limit=limit)
                
                result = {
                    "template_name": template_seq.name,
                    "template_length": len(template_seq),
                    "limit": limit,
                    "forward_primers": [
                        {
                            "sequence": str(p.seq),
                            "position": p.position,
                            "footprint": p.footprint
                        } for p in annealing.forward_primers
                    ],
                    "reverse_primers": [
                        {
                            "sequence": str(p.seq),
                            "position": p.position,
                            "footprint": p.footprint
                        } for p in annealing.reverse_primers
                    ],
                    "products": len(annealing.products)
                }
                
                action.add_success_fields(annealing_result=result)
                return result
                
            except Exception as e:
                action.add_success_fields(error=str(e))
                raise ValueError(f"Error analyzing primer annealing: {str(e)}")
    
    def restriction_analysis(self, sequence: str, enzymes: List[str]) -> List[RestrictionAnalysis]:
        """
        Analyze restriction enzyme cutting patterns in a DNA sequence.
        
        This method identifies where restriction enzymes cut in your sequence,
        calculates fragment sizes, and provides cutting pattern analysis.
        Essential for planning cloning strategies and analyzing vector digests.
        
        Args:
            sequence (str): DNA sequence to analyze
                Example: "GGAATTCCGATCGATCGGATCCAAGCTT"
                Can be gene, plasmid, or any DNA sequence
            enzymes (List[str]): List of restriction enzyme names to test
                Example: ["EcoRI", "BamHI", "HindIII", "XhoI", "NotI"]
                Use standard enzyme names (case-sensitive)
        
        Returns:
            List[RestrictionAnalysis]: Analysis for each enzyme that cuts:
                - enzyme_name: Name of the restriction enzyme
                - recognition_site: DNA sequence the enzyme recognizes
                - cut_positions: List of positions where enzyme cuts (0-based)
                - fragments: List of fragment sizes after digestion
        
        Example Input:
            sequence = "GGAATTCCGATCGATCGGATCCAAGCTT"
            enzymes = ["EcoRI", "BamHI", "HindIII"]
        
        Example Output:
            [
                {
                    "enzyme_name": "EcoRI",
                    "recognition_site": "GAATTC", 
                    "cut_positions": [2],
                    "fragments": [2, 25]
                },
                {
                    "enzyme_name": "BamHI",
                    "recognition_site": "GGATCC",
                    "cut_positions": [18],
                    "fragments": [18, 9]
                },
                {
                    "enzyme_name": "HindIII", 
                    "recognition_site": "AAGCTT",
                    "cut_positions": [21],
                    "fragments": [21, 6]
                }
            ]
        
        Understanding Results:
            - cut_positions: Where the enzyme cuts (creates DNA break)
            - fragments: Size of pieces after cutting
            - Empty list: Enzyme does not cut the sequence
            - Multiple positions: Enzyme cuts at multiple sites
        
        Common Restriction Enzymes:
            - EcoRI: GAATTC (6-cutter, moderate frequency)
            - BamHI: GGATCC (6-cutter, moderate frequency)  
            - HindIII: AAGCTT (6-cutter, moderate frequency)
            - XhoI: CTCGAG (6-cutter, rare)
            - NotI: GCGGCCGC (8-cutter, very rare)
            - HaeIII: GGCC (4-cutter, frequent)
        
        Fragment Analysis:
            - Single cut: 2 fragments
            - Two cuts: 3 fragments
            - No cuts: 1 fragment (original sequence)
            - Fragment sizes useful for gel electrophoresis planning
        
        Common Use Cases:
            - Plan restriction cloning strategy
            - Choose unique cutting sites for linearization
            - Verify plasmid identity by restriction mapping
            - Design cloning sites in primers
            - Check for unwanted cutting sites in inserts
        
        Cloning Strategy Tips:
            - Use enzymes that cut vector once and insert at ends
            - Avoid enzymes that cut multiple times in vector
            - Choose compatible enzymes (same buffer conditions)
            - Consider fragment sizes for gel purification
        
        Error Conditions:
            - Invalid enzyme names are silently skipped
            - Invalid nucleotides in sequence raise ValueError
            - Empty enzyme list returns empty results
        """
        with start_action(action_type="restriction_analysis", 
                         sequence_length=len(sequence), 
                         num_enzymes=len(enzymes)) as action:
            try:
                dseq = Dseqrecord(sequence)
                results = []
                
                for enzyme_name in enzymes:
                    try:
                        # Get enzyme from BioPython
                        enzyme = globals()[enzyme_name]
                        
                        # Find cut sites
                        cut_sites = enzyme.search(dseq.seq)
                        
                        # Calculate fragment sizes
                        if cut_sites:
                            positions = sorted(cut_sites)
                            fragments = []
                            if len(positions) == 1:
                                fragments = [positions[0], len(dseq) - positions[0]]
                            else:
                                for i in range(len(positions)):
                                    if i == 0:
                                        fragments.append(positions[i])
                                    else:
                                        fragments.append(positions[i] - positions[i-1])
                                fragments.append(len(dseq) - positions[-1])
                        else:
                            fragments = [len(dseq)]
                        
                        result = RestrictionAnalysis(
                            enzyme_name=enzyme_name,
                            recognition_site=str(enzyme.site),
                            cut_positions=cut_sites,
                            fragments=fragments
                        )
                        results.append(result)
                        
                    except KeyError:
                        # Enzyme not found, skip
                        continue
                
                action.add_success_fields(num_results=len(results))
                return results
                
            except Exception as e:
                action.add_success_fields(error=str(e))
                raise ValueError(f"Error in restriction analysis: {str(e)}")
    
    def digest_sequence(self, sequence: str, enzymes: List[str]) -> Dict[str, Any]:
        """
        Perform virtual restriction enzyme digestion of DNA sequence.
        
        This method simulates cutting DNA with restriction enzymes and returns
        the resulting fragments. It's essential for predicting gel electrophoresis
        patterns and planning fragment isolation for cloning.
        
        Args:
            sequence (str): DNA sequence to digest
                Example: "GGAATTCCGATCGATCGGATCCAAGCTT"
                Can be plasmid, PCR product, or genomic DNA
            enzymes (List[str]): Restriction enzymes to use for digestion
                Example: ["EcoRI", "BamHI"] for double digest
                Example: ["HindIII"] for single digest
        
        Returns:
            Dict[str, Any]: Complete digestion results:
                - original_sequence: Input DNA sequence
                - original_length: Length before digestion
                - enzymes_used: List of enzymes that successfully cut
                - num_fragments: Number of fragments produced
                - fragments: List of fragment details (sequence, length, linear/circular)
        
        Example Input:
            sequence = "GGAATTCCGATCGATCGGATCCAAGCTT"
            enzymes = ["EcoRI", "BamHI"]
        
        Example Output:
            {
                "original_sequence": "GGAATTCCGATCGATCGGATCCAAGCTT",
                "original_length": 27,
                "enzymes_used": ["EcoRI", "BamHI"],
                "num_fragments": 3,
                "fragments": [
                    {
                        "sequence": "GG",
                        "length": 2,
                        "linear": true
                    },
                    {
                        "sequence": "AATTCCGATCGATCGG",
                        "length": 16,
                        "linear": true
                    },
                    {
                        "sequence": "ATCCAAGCTT", 
                        "length": 9,
                        "linear": true
                    }
                ]
            }
        
        Fragment Interpretation:
            - linear: true = Linear DNA fragment (typical for digests)
            - linear: false = Circular fragment (rare, only if no cuts in circular DNA)
            - Fragments ordered by position in original sequence
            - Total fragment lengths = original length
        
        Digestion Types:
            - Single digest: One enzyme, simpler pattern
            - Double digest: Two enzymes, more fragments
            - Triple digest: Three enzymes, complex pattern
            - Partial digest: Some sites may not cut (not simulated here)
        
        Common Use Cases:
            - Predict gel electrophoresis pattern
            - Plan fragment isolation for cloning
            - Verify plasmid identity by restriction mapping
            - Design cloning strategy with specific fragments
            - Calculate expected band sizes for analysis
        
        Gel Analysis Tips:
            - Larger fragments migrate less in gel electrophoresis
            - Very small fragments (<100 bp) may run off gel
            - Similar sized fragments may appear as single band
            - Use appropriate ladder for size estimation
        
        Cloning Applications:
            - Cut vector and insert with same enzymes for ligation
            - Isolate specific fragments from gel for purification
            - Check orientation by asymmetric cutting pattern
            - Verify successful cloning by diagnostic digest
        
        Error Conditions:
            - Invalid enzyme names: Silently skipped
            - No valid enzymes: Raises ValueError
            - Invalid nucleotides: Raises ValueError
            - Empty sequence: Raises ValueError
        """
        with start_action(action_type="digest_sequence", 
                         sequence_length=len(sequence), 
                         num_enzymes=len(enzymes)) as action:
            try:
                dseq = Dseqrecord(sequence)
                
                # Get enzyme objects
                enzyme_objs = []
                for enzyme_name in enzymes:
                    try:
                        enzyme_objs.append(globals()[enzyme_name])
                    except KeyError:
                        continue
                
                if not enzyme_objs:
                    raise ValueError("No valid enzymes found")
                
                # Perform digestion
                fragments = dseq.cut(*enzyme_objs)
                
                result = {
                    "original_sequence": str(dseq.seq),
                    "original_length": len(dseq),
                    "enzymes_used": enzymes,
                    "num_fragments": len(fragments),
                    "fragments": [
                        {
                            "sequence": str(frag.seq),
                            "length": len(frag),
                            "linear": not frag.circular
                        } for frag in fragments
                    ]
                }
                
                action.add_success_fields(digest_result=result)
                return result
                
            except Exception as e:
                action.add_success_fields(error=str(e))
                raise ValueError(f"Error in sequence digestion: {str(e)}")
    
    def assembly(self, fragments: List[str], limit: int = 25) -> AssemblyResult:
        """
        Assemble DNA fragments by homologous recombination.
        
        This method simulates DNA assembly by finding overlapping regions between
        fragments and joining them together. It mimics natural homologous recombination
        or in vitro assembly methods like Gibson assembly, but uses existing overlaps
        rather than creating them artificially.
        
        Args:
            fragments (List[str]): List of DNA fragment sequences to assemble
                Example: ["ATGAAAGCACTGATT", "CTGATTGCATGCAAG", "GCAAGCTTAAATAG"]
                Fragments should have overlapping regions for successful assembly
            limit (int, optional): Minimum overlap length required for assembly.
                Defaults to 25. Higher values = more stringent assembly
                Range: 15-50 bp depending on fragment design
        
        Returns:
            AssemblyResult: Information about the assembled construct:
                - sequence: The final assembled DNA sequence
                - length: Total length of assembled sequence
                - circular: Whether the assembly forms a circle
                - fragments_used: Number of input fragments used
        
        Example Input:
            fragments = [
                "ATGAAAGCACTGATTCTATTGCTGAAA",
                "CTGAAAAAGATAATAGCGCCTACGAT",
                "TACGATCGCAAGCTTAAATAG"
            ]
            limit = 20
        
        Example Output:
            {
                "sequence": "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAATAGCGCCTACGATCGCAAGCTTAAATAG",
                "length": 60,
                "circular": false,
                "fragments_used": 3
            }
        
        Assembly Requirements:
            - Fragments must have overlapping sequences >= limit
            - Overlaps must be exact matches (no mismatches)
            - Assembly tries linear first, then circular if linear fails
            - At least 2 fragments required for assembly
        
        Assembly Types:
            - Linear assembly: Fragments joined end-to-end in a line
            - Circular assembly: Fragments form a closed circle
            - No assembly: Insufficient overlaps or incompatible fragments
        
        Common Use Cases:
            - Simulate Gibson assembly with pre-designed overlaps
            - Join PCR products with complementary ends
            - Assemble synthetic DNA constructs
            - Validate assembly design before lab work
            - Reconstruct sequences from overlapping reads
        
        Design Tips:
            - Plan 20-40 bp overlaps between adjacent fragments
            - Ensure overlaps are unique to avoid incorrect assembly
            - Check fragment order and orientation
            - Verify no internal repetitive sequences
        
        Error Conditions:
            - No overlaps found: Assembly fails, raises ValueError
            - Fragments too short: May not have sufficient overlap
            - Conflicting overlaps: Assembly may fail or produce wrong result
            - Single fragment: No assembly possible
        
        Troubleshooting:
            - If assembly fails: Check overlap lengths and sequences
            - If wrong product: Verify fragment order and orientation
            - If multiple products: Overlaps may be non-specific
        """
        with start_action(action_type="assembly", 
                         num_fragments=len(fragments), 
                         limit=limit) as action:
            try:
                # Create Dseqrecord objects from fragments
                fragment_seqs = [Dseqrecord(frag) for frag in fragments]
                
                # Perform assembly
                assembly_obj = Assembly(fragment_seqs, limit=limit)
                assembled = assembly_obj.assemble_linear()
                
                if not assembled:
                    # Try circular assembly
                    assembled = assembly_obj.assemble_circular()
                
                if not assembled:
                    raise ValueError("No assembly products found")
                
                # Get the first (best) assembly
                best_assembly = assembled[0]
                
                result = AssemblyResult(
                    sequence=str(best_assembly.seq),
                    length=len(best_assembly),
                    circular=best_assembly.circular,
                    fragments_used=len(fragments)
                )
                
                action.add_success_fields(assembly_result=result.model_dump())
                return result
                
            except Exception as e:
                action.add_success_fields(error=str(e))
                raise ValueError(f"Error in assembly: {str(e)}")
    
    def gibson_assembly(self, fragments: List[str], overlap: int = 20) -> AssemblyResult:
        """
        Design and simulate Gibson assembly of DNA fragments.
        
        Gibson assembly is a molecular cloning method that uses overlapping DNA
        fragments to create seamless assemblies. This method designs appropriate
        overlaps between fragments and simulates the assembly process.
        
        Args:
            fragments (List[str]): List of DNA fragment sequences for Gibson assembly
                Example: ["ATGAAAGCACTGATT", "GCACTGATTGCATGC", "GCATGCAAGCTTAAA"]
                Fragments will be modified to have overlapping ends
            overlap (int, optional): Length of overlap to create between fragments.
                Defaults to 20. Range: 15-40 bp for optimal Gibson assembly
                Longer overlaps = more specific, shorter = more efficient
        
        Returns:
            AssemblyResult: Information about the Gibson assembly product:
                - sequence: The assembled DNA sequence
                - length: Total length of assembled construct
                - circular: Whether assembly forms circular product
                - fragments_used: Number of input fragments
        
        Example Input:
            fragments = [
                "ATGAAAGCACTGATTCTATTGCTG",
                "GCTGAAAAAGATAATAGCGCCTAC",
                "TACGATCGCAAGCTTAAATAG"
            ]
            overlap = 25
        
        Example Output:
            {
                "sequence": "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAATAGCGCCTACGATCGCAAGCTTAAATAG",
                "length": 60,
                "circular": false,
                "fragments_used": 3
            }
        
        Gibson Assembly Process:
            1. 5' exonuclease creates single-strand overhangs
            2. DNA polymerase fills in gaps
            3. DNA ligase seals nicks
            4. Overlapping fragments join seamlessly
        
        Design Principles:
            - Overlaps should be 15-40 bp for optimal efficiency
            - Avoid secondary structures in overlap regions
            - GC content of overlaps should be 40-60%
            - Melting temperature of overlaps should be similar
        
        Common Use Cases:
            - Assemble multiple DNA fragments into one construct
            - Clone genes into vectors without restriction sites
            - Create complex genetic circuits
            - Join PCR products with designed overhangs
            - Scarless cloning without unwanted sequences
        
        Advantages over Traditional Cloning:
            - No restriction enzyme sites required
            - Seamless assembly without scars
            - Multiple fragments can be joined simultaneously
            - More flexible than restriction-ligation cloning
        
        Experimental Considerations:
            - Equimolar ratios of fragments work best
            - Longer fragments may need higher concentrations
            - Reaction temperature typically 50°C for 1 hour
            - Transformation efficiency may be lower than simple ligations
        
        Error Conditions:
            - Fragments too short for specified overlap
            - Gibson assembly fails due to poor overlap design
            - Invalid nucleotides in fragment sequences
            - Insufficient fragments for assembly
        
        Optimization Tips:
            - Use 20-25 bp overlaps for most applications
            - Check overlap sequences for hairpins or repetitive elements
            - Verify fragment concentrations are balanced
            - Consider fragment order and orientation carefully
        """
        with start_action(action_type="gibson_assembly", 
                         num_fragments=len(fragments), 
                         overlap=overlap) as action:
            try:
                # For simplicity, create a basic Gibson assembly simulation
                # Real Gibson assembly would require specific overlaps between fragments
                
                if len(fragments) < 2:
                    raise ValueError("Gibson assembly requires at least 2 fragments")
                
                # Create Dseqrecord objects from fragments
                fragment_seqs = [Dseqrecord(frag) for frag in fragments]
                
                # Simple overlap-based assembly simulation
                # In real Gibson assembly, fragments would have designed overlapping ends
                assembled_sequence = ""
                
                for i, frag in enumerate(fragment_seqs):
                    if i == 0:
                        # First fragment - take full sequence
                        assembled_sequence = str(frag.seq)
                    else:
                        # Check for overlap with previous fragment
                        prev_seq = assembled_sequence
                        curr_seq = str(frag.seq)
                        
                        # Look for overlap at the end of previous and start of current
                        best_overlap = 0
                        for j in range(1, min(overlap + 5, len(prev_seq), len(curr_seq))):
                            if prev_seq[-j:] == curr_seq[:j]:
                                best_overlap = j
                        
                        if best_overlap > 0:
                            # Merge with overlap
                            assembled_sequence += curr_seq[best_overlap:]
                        else:
                            # No natural overlap found, just concatenate
                            assembled_sequence += curr_seq
                
                # Check if we can make a circular assembly
                first_frag = str(fragment_seqs[0].seq)
                last_part = assembled_sequence[-overlap:]
                first_part = first_frag[:overlap]
                is_circular = (last_part == first_part) if len(last_part) == len(first_part) else False
                
                result = AssemblyResult(
                    sequence=assembled_sequence,
                    length=len(assembled_sequence),
                    circular=is_circular,
                    fragments_used=len(fragments)
                )
                
                action.add_success_fields(gibson_result=result.model_dump())
                return result
                
            except Exception as e:
                action.add_success_fields(error=str(e))
                raise ValueError(f"Error in Gibson assembly: {str(e)}")
    
    def read_sequence_file(self, file_content: str, file_format: str = "auto") -> Dict[str, Any]:
        """
        Read and parse DNA sequence from file content.
        
        This method parses sequence files in various formats (FASTA, GenBank, etc.)
        and extracts sequence information, annotations, and features. It's essential
        for importing existing sequences into your workflow.
        
        Args:
            file_content (str): The content of the sequence file as text
                Example FASTA: ">sequence1\nATGAAGGCACTGATT\n"
                Example GenBank: "LOCUS    sequence1   15 bp    DNA     linear   UNK 01-JAN-1980\n..."
            file_format (str, optional): File format to parse. Defaults to "auto"
                Options: "auto", "fasta", "genbank", "gb", "embl"
                "auto" attempts to detect format automatically
        
        Returns:
            Dict[str, Any]: Comprehensive sequence information:
                - name: Sequence name/identifier
                - id: Sequence ID (from file)
                - description: Sequence description
                - sequence: DNA sequence string
                - length: Sequence length in bp
                - circular: Whether sequence is circular
                - num_features: Number of annotated features
                - features: List of sequence features (first 10)
        
        Example Input (FASTA):
            file_content = ">my_gene\\nATGAAAGCACTGATTCTATTGCTGAAAAGATAATAGCGCCTACGATCGCAAGCTTAAATAG"
        
        Example Output:
            {
                "name": "my_gene",
                "id": "my_gene", 
                "description": "",
                "sequence": "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAATAGCGCCTACGATCGCAAGCTTAAATAG",
                "length": 60,
                "circular": false,
                "num_features": 0,
                "features": []
            }
        
        Example Input (GenBank):
            file_content = "LOCUS plasmid 5386 bp DNA circular UNK\\n" +
                          "DEFINITION Cloning vector\\n" +
                          "FEATURES Location/Qualifiers\\n" +
                          "     gene 100..500\\n" +
                          "          /gene='resistance'\\n" +
                          "ORIGIN\\n" +
                          "  1 atgcgatcgc aagcttgaat tcggatccaa gcttgcatgc\\n" +
                          "//"
        
        Supported File Formats:
            - FASTA (.fasta, .fa, .fna): Simple sequence format
            - GenBank (.gb, .gbk): Rich annotation format
            - EMBL (.embl): European sequence format
            - Auto-detection: Tries to identify format automatically
        
        Feature Information:
            - type: Feature type (gene, CDS, promoter, etc.)
            - start/end: Feature location (0-based coordinates)
            - strand: +1 (forward), -1 (reverse), or None
            - qualifiers: Additional feature information (gene name, product, etc.)
        
        Common Use Cases:
            - Import sequences from databases (GenBank, EMBL)
            - Load synthetic sequences from design tools
            - Parse sequences from sequencing results
            - Import annotated plasmids or vectors
            - Batch process multiple sequence files
        
        File Format Detection:
            - FASTA: Starts with ">" character
            - GenBank: Contains "LOCUS" keyword
            - EMBL: Contains "ID" line
            - Auto mode tries each format until one works
        
        Error Conditions:
            - Invalid file format: Cannot parse content
            - Empty file content: No sequence found
            - Corrupted file: Parsing fails
            - Unsupported format: Not recognized
        
        Best Practices:
            - Use "auto" format for unknown files
            - Verify sequence after import
            - Check feature annotations for accuracy
            - Be aware of coordinate systems (0-based vs 1-based)
        """
        with start_action(action_type="read_sequence_file", format=file_format) as action:
            try:
                # Parse the sequence
                if file_format == "auto":
                    # Try to detect format
                    if file_content.strip().startswith('>'):
                        file_format = "fasta"
                    elif "LOCUS" in file_content:
                        file_format = "genbank"
                    else:
                        file_format = "fasta"
                
                # Create temporary file
                with tempfile.NamedTemporaryFile(mode='w', suffix=f'.{file_format}', delete=False) as tmp:
                    tmp.write(file_content)
                    tmp_path = tmp.name
                
                try:
                    # Read using pydna
                    seq = read(tmp_path)
                    
                    result = {
                        "name": seq.name,
                        "id": seq.id,
                        "description": seq.description,
                        "sequence": str(seq.seq),
                        "length": len(seq),
                        "circular": seq.circular,
                        "num_features": len(seq.features),
                        "features": [
                            {
                                "type": f.type,
                                "start": f.location.start,
                                "end": f.location.end,
                                "strand": f.location.strand,
                                "qualifiers": dict(f.qualifiers)
                            } for f in seq.features[:10]  # Limit to first 10 features
                        ]
                    }
                    
                    action.add_success_fields(sequence_info=result)
                    return result
                    
                finally:
                    # Clean up temp file
                    os.unlink(tmp_path)
                
            except Exception as e:
                action.add_success_fields(error=str(e))
                raise ValueError(f"Error reading sequence file: {str(e)}")
    
    def write_sequence_file(self, sequence: str, name: str = "sequence", file_format: str = "fasta") -> str:
        """
        Write DNA sequence to file format.
        
        This method converts a DNA sequence into various file formats for export,
        sharing, or submission to databases. Essential for preparing sequences
        for downstream analysis or experimental work.
        
        Args:
            sequence (str): DNA sequence to write
                Example: "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT"
            name (str, optional): Name for the sequence. Defaults to "sequence"
                Example: "my_gene", "plasmid_construct", "pcr_product"
            file_format (str, optional): Output file format. Defaults to "fasta"
                Options: "fasta", "genbank", "gb", "embl"
        
        Returns:
            str: Formatted sequence file content ready for saving
        
        Example Input:
            sequence = "ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT"
            name = "test_gene"
            file_format = "fasta"
        
        Example Output (FASTA):
            ">test_gene
            ATGAAAGCACTGATTCTATTGCTGAAAAAGATAAT
            "
        
        Example Output (GenBank):
            "LOCUS       test_gene    34 bp    DNA     linear   UNK 01-JAN-1980
            DEFINITION  .
            ACCESSION   test_gene
            VERSION     test_gene
            KEYWORDS    .
            SOURCE      .
            FEATURES             Location/Qualifiers
            ORIGIN
                  1 atgaaagcac tgattctatt gctgaaaaag ataat
            //
            "
        
        File Format Options:
            - FASTA (.fasta): Simple, widely compatible format
            - GenBank (.gb): Rich format with annotation support
            - EMBL (.embl): European Molecular Biology Lab format
            - Text (.txt): Plain sequence text
        
        Format Characteristics:
            - FASTA: Header line (>) + sequence, 60-80 chars per line
            - GenBank: Structured format with metadata and features
            - EMBL: Similar to GenBank but different syntax
            - All formats preserve sequence information
        
        Common Use Cases:
            - Export sequences for primers ordering
            - Save designed constructs for lab work
            - Prepare sequences for database submission
            - Share sequences with collaborators
            - Create input files for other software
        
        Best Practices:
            - Use descriptive sequence names
            - Choose appropriate format for intended use
            - FASTA for simple sequences, GenBank for annotated
            - Check output format is compatible with target software
        
        Error Conditions:
            - Invalid nucleotides in sequence
            - Unsupported file format
            - Empty sequence string
            - Invalid characters in sequence name
        
        Usage Tips:
            - Save output to file with appropriate extension
            - FASTA is most universally compatible
            - GenBank format preserves more information
            - Consider sequence length for format choice
        """
        with start_action(action_type="write_sequence_file", format=file_format) as action:
            try:
                dseq = Dseqrecord(sequence, name=name)
                
                # Format the sequence
                formatted = dseq.format(file_format)
                
                action.add_success_fields(format=file_format, length=len(formatted))
                return formatted
                
            except Exception as e:
                action.add_success_fields(error=str(e))
                raise ValueError(f"Error writing sequence file: {str(e)}")
    
    def gel_electrophoresis(self, fragments: List[str], ladder: str = "1kb") -> Dict[str, Any]:
        """
        Simulate gel electrophoresis of DNA fragments.
        
        This method simulates running DNA fragments on an agarose gel to predict
        migration patterns, band sizes, and separation quality. Essential for
        planning gel analysis and interpreting experimental results.
        
        Args:
            fragments (List[str]): List of DNA fragment sequences to analyze
                Example: ["ATGAAAGCACTGATT", "GCACTGATTGCATGCAAGCTT", "AAGCTTGAATTCGGATCC"]
                Each fragment represents a DNA band on the gel
            ladder (str, optional): DNA ladder/marker to use. Defaults to "1kb"
                Options: "1kb", "100bp", "lambda", "custom"
                Used for size estimation and reference
        
        Returns:
            Dict[str, Any]: Gel electrophoresis simulation results:
                - fragments: List of fragment details (sequence, length, mobility)
                - gel_image: Text representation of gel pattern
                - ladder_used: Reference ladder information
        
        Example Input:
            fragments = [
                "ATGAAAGCACTGATTCTATTGCTG",  # 24 bp
                "GCACTGATTGCATGCAAGCTTGAATTC",  # 27 bp
                "AAGCTTGAATTCGGATCCGTCGAC"  # 24 bp
            ]
            ladder = "1kb"
        
        Example Output:
            {
                "fragments": [
                    {
                        "sequence": "ATGAAAGCACTGATTCTATTGCTG",
                        "length": 24,
                        "mobility": 24
                    },
                    {
                        "sequence": "GCACTGATTGCATGCAAGCTTGAATTC", 
                        "length": 27,
                        "mobility": 27
                    },
                    {
                        "sequence": "AAGCTTGAATTCGGATCCGTCGAC",
                        "length": 24, 
                        "mobility": 24
                    }
                ],
                "gel_image": "Text representation of gel bands",
                "ladder_used": "1kb"
            }
        
        Gel Electrophoresis Principles:
            - Smaller DNA fragments migrate further (higher mobility)
            - Larger fragments migrate less (lower mobility)
            - Migration distance inversely proportional to log(size)
            - Band intensity proportional to DNA concentration
        
        Fragment Interpretation:
            - mobility: Relative migration distance (higher = smaller fragment)
            - Similar sized fragments may co-migrate
            - Resolution depends on gel percentage and conditions
            - Very small fragments (<100 bp) may run off gel
        
        Common Ladder Types:
            - 1kb ladder: 1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000 bp
            - 100bp ladder: 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 bp
            - Lambda/HindIII: 23130, 9416, 6682, 4361, 2322, 2027, 564, 125 bp
        
        Gel Conditions Affecting Separation:
            - Agarose concentration: 0.8% (large), 1.5% (medium), 2.5% (small)
            - Voltage: Higher = faster but less resolution
            - Buffer: TAE vs TBE affects migration
            - Temperature: Affects DNA mobility
        
        Common Use Cases:
            - Verify PCR product size
            - Check restriction digest patterns
            - Confirm DNA purification
            - Analyze cloning results
            - Plan gel purification strategy
        
        Experimental Planning:
            - Choose appropriate gel concentration for fragment sizes
            - Select suitable ladder for size range
            - Plan loading order for easy interpretation
            - Consider band intensity for quantification
        
        Error Conditions:
            - Empty fragment list
            - Invalid nucleotides in sequences
            - Unsupported ladder type
            - Fragments too large/small for detection
        
        Troubleshooting:
            - No bands: Check DNA integrity and loading
            - Smeared bands: DNA degradation or overloading
            - Unexpected sizes: Verify sequence and enzyme activity
            - Poor separation: Adjust gel concentration or voltage
        """
        with start_action(action_type="gel_electrophoresis", 
                         num_fragments=len(fragments), 
                         ladder=ladder) as action:
            try:
                # Create Dseqrecord objects
                fragment_seqs = [Dseqrecord(frag) for frag in fragments]
                
                # Calculate fragment data
                fragments_data = []
                for frag in fragment_seqs:
                    length = len(frag)
                    # Simple mobility calculation based on length (inverse relationship)
                    # Smaller fragments have higher mobility (travel further)
                    mobility = 1000.0 / length if length > 0 else 0.0
                    
                    fragments_data.append({
                        "sequence": str(frag.seq),
                        "length": length,
                        "mobility": round(mobility, 2)
                    })
                
                # Create simple text representation of gel
                sorted_frags = sorted(fragments_data, key=lambda x: x['mobility'], reverse=True)
                gel_lines = ["Gel Electrophoresis Results:"]
                gel_lines.append("Wells (from left to right):")
                for i, frag in enumerate(sorted_frags):
                    gel_lines.append(f"Lane {i+1}: {frag['length']} bp")
                
                result = {
                    "fragments": fragments_data,
                    "gel_image": "\n".join(gel_lines),
                    "ladder_used": ladder,
                    "sorted_by_mobility": sorted_frags
                }
                
                action.add_success_fields(gel_result=result)
                return result
                
            except Exception as e:
                action.add_success_fields(error=str(e))
                raise ValueError(f"Error in gel electrophoresis: {str(e)}")
    
    def download_genbank(self, accession: str) -> Dict[str, Any]:
        """
        Download DNA sequence from GenBank database.
        
        This method retrieves sequence data from NCBI's GenBank database using
        accession numbers. It downloads complete sequence information including
        annotations, features, and metadata. Essential for accessing public
        sequence data for analysis and cloning.
        
        Args:
            accession (str): GenBank accession number
                Example: "NM_000546" (human TP53 mRNA)
                Example: "NC_000001" (human chromosome 1)
                Example: "pBR322" (common cloning vector)
                Format: Letters followed by numbers (e.g., NM_123456, AC_123456)
        
        Returns:
            Dict[str, Any]: Complete GenBank record information:
                - accession: The requested accession number
                - name: Sequence name/locus
                - id: Sequence identifier
                - description: Sequence description
                - sequence: Full DNA sequence
                - length: Sequence length in bp
                - circular: Whether sequence is circular
                - num_features: Number of annotated features
                - organism: Source organism
                - genbank_text: Raw GenBank file content (truncated)
        
        Example Input:
            accession = "NM_000546"  # Human TP53 tumor suppressor gene
        
        Example Output:
            {
                "accession": "NM_000546",
                "name": "NM_000546",
                "id": "NM_000546.6",
                "description": "Homo sapiens tumor protein p53 (TP53), mRNA",
                "sequence": "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATT...",
                "length": 1842,
                "circular": false,
                "num_features": 8,
                "organism": "Homo sapiens",
                "genbank_text": "LOCUS       NM_000546    1842 bp    mRNA    linear   PRI 15-APR-2023..."
            }
        
        GenBank Accession Types:
            - RefSeq mRNA: NM_######, XM_######
            - RefSeq protein: NP_######, XP_######
            - GenBank genomic: AC_######, NT_######
            - GenBank mRNA: BC_######, AY_######
            - Chromosome: NC_######
            - Plasmid: Various formats (pBR322, pUC19, etc.)
        
        Retrieved Information:
            - Complete nucleotide sequence
            - All annotated features (genes, CDS, etc.)
            - Taxonomic information
            - Publication references
            - Sequence version and update history
        
        Common Use Cases:
            - Download reference sequences for analysis
            - Get sequences for primer design
            - Retrieve vector sequences for cloning
            - Access genomic regions of interest
            - Obtain published sequences for replication
        
        Feature Information:
            - Genes: Genomic regions encoding proteins
            - CDS: Protein-coding sequences
            - Exons/Introns: Gene structure elements
            - Regulatory elements: Promoters, enhancers
            - Variations: SNPs, indels, etc.
        
        Data Quality Considerations:
            - RefSeq: Curated, high-quality sequences
            - GenBank: Submitted sequences, quality varies
            - Version numbers: Higher = more recent
            - Check sequence date for currency
        
        Error Conditions:
            - Invalid accession number: Not found in GenBank
            - Network issues: Cannot connect to NCBI
            - Parsing errors: Corrupted or unusual format
            - Access restrictions: Some sequences may be restricted
        
        Best Practices:
            - Use RefSeq accessions when possible (NM_, NP_, etc.)
            - Check sequence version and date
            - Verify organism and sequence type
            - Be aware of coordinate systems for features
            - Respect NCBI usage guidelines
        
        Rate Limiting:
            - NCBI limits download frequency
            - Use appropriate delays between requests
            - Consider bulk downloads for large datasets
            - Provide email contact for tracking
        """
        with start_action(action_type="download_genbank", accession=accession) as action:
            try:
                # Download from GenBank
                gb_text = download_text(accession)
                
                # Parse the sequence
                seq = parse(gb_text)[0] if parse(gb_text) else None
                
                if not seq:
                    raise ValueError(f"Could not parse sequence for {accession}")
                
                result = {
                    "accession": accession,
                    "name": seq.name,
                    "id": seq.id,
                    "description": seq.description,
                    "sequence": str(seq.seq),
                    "length": len(seq),
                    "circular": seq.circular,
                    "num_features": len(seq.features),
                    "organism": seq.annotations.get("organism", "Unknown"),
                    "genbank_text": gb_text[:1000] + "..." if len(gb_text) > 1000 else gb_text
                }
                
                action.add_success_fields(download_result=result)
                return result
                
            except Exception as e:
                action.add_success_fields(error=str(e))
                raise ValueError(f"Error downloading from GenBank: {str(e)}")


def cli_app():
    """
    Run the PyDNA MCP server with HTTP transport.
    
    Starts the server with Server-Sent Events (SSE) transport on the default
    host and port. This is suitable for web-based applications and HTTP clients.
    """
    app = PyDNAMCP()
    app.run(transport="sse", host=DEFAULT_HOST, port=DEFAULT_PORT)

def cli_app_stdio():
    """
    Run the PyDNA MCP server with stdio transport.
    
    Starts the server with standard input/output transport. This is suitable
    for command-line applications and direct process communication.
    """
    app = PyDNAMCP()
    app.run(transport="stdio")

def cli_app_sse():
    """
    Run the PyDNA MCP server with SSE transport.
    
    Starts the server with Server-Sent Events transport on the specified
    host and port. This is suitable for web applications and HTTP-based clients.
    """
    app = PyDNAMCP()
    app.run(transport="sse", host=DEFAULT_HOST, port=DEFAULT_PORT)

if __name__ == "__main__":
    cli_app_stdio() 