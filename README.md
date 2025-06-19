# pydna-mcp

[![Tests](https://github.com/antonkulaga/pydna-mcp/actions/workflows/test.yml/badge.svg)](https://github.com/antonkulaga/pydna-mcp/actions/workflows/test.yml)
[![PyPI version](https://badge.fury.io/py/pydna-mcp.svg)](https://badge.fury.io/py/pydna-mcp)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

MCP (Model Context Protocol) server for PyDNA - comprehensive DNA simulation and molecular cloning tools

This server implements the Model Context Protocol (MCP) for PyDNA, providing a standardized interface for DNA simulation and molecular cloning operations. MCP enables AI assistants and agents to perform sophisticated molecular biology tasks through structured interfaces, making computational DNA manipulation as natural as describing your experimental workflow.

The server wraps the powerful PyDNA library to offer:

- **PCR Simulation**: Design primers, simulate amplification, analyze products
- **Primer Design**: Automated primer design with customizable parameters  
- **Restriction Analysis**: Enzyme mapping, digestion simulation, fragment analysis
- **DNA Assembly**: Gibson assembly, homologous recombination, fragment joining
- **Sequence Manipulation**: Creation, analysis, reverse complement, format conversion
- **Gel Electrophoresis**: Simulate gel runs for fragment size analysis
- **GenBank Integration**: Download sequences, read/write various file formats
- **Molecular Calculations**: Melting temperatures, GC content, molecular weights

If you want to understand more about what the Model Context Protocol is and how to use it more efficiently, you can take the [DeepLearning AI Course](https://www.deeplearning.ai/short-courses/mcp-build-rich-context-ai-apps-with-anthropic/) or search for MCP videos on YouTube.

## 🏆 Part of Holy Bio MCP Framework  

This MCP server is part of the **[Holy Bio MCP](https://github.com/longevity-genie/holy-bio-mcp)** project - a unified framework for bioinformatics research that **won the Bio x AI Hackathon 2025** and continues to be actively developed and extended after the victory.

The Holy Bio MCP framework brings together multiple specialized MCP servers into a cohesive ecosystem for advanced biological research:

- **[gget-mcp](https://github.com/longevity-genie/gget-mcp)** - Genomics & sequence analysis toolkit
- **[opengenes-mcp](https://github.com/longevity-genie/opengenes-mcp)** - Aging & longevity genetics
- **[synergy-age-mcp](https://github.com/longevity-genie/synergy-age-mcp)** - Synergistic genetic interactions in longevity
- **[biothings-mcp](https://github.com/longevity-genie/biothings-mcp)** - Foundational biological data from BioThings.io
- **[pharmacology-mcp](https://github.com/antonkulaga/pharmacology-mcp)** - Drug, target & ligand data
- **[pydna-mcp](https://github.com/antonkulaga/pydna-mcp)** - DNA simulation & molecular cloning (this server)

Together, these servers provide 50+ specialized bioinformatics functions that can work seamlessly together in AI-driven research workflows. Learn more about the complete framework at [github.com/longevity-genie/holy-bio-mcp](https://github.com/longevity-genie/holy-bio-mcp).

## Usage Example

Here's how PyDNA MCP server enables natural language molecular cloning:

```
User: "I want to amplify the GFP gene from pGLO plasmid using primers with EcoRI sites for cloning into pUC19"

AI Assistant: I'll help you design that cloning strategy:
1. First, let me design primers with EcoRI sites for GFP
2. Simulate the PCR amplification  
3. Analyze the restriction sites in both your PCR product and pUC19
4. Show you the expected fragments after digestion

[Uses pydna-mcp tools to execute each step]
```

*The server translates complex molecular biology workflows into simple tool calls, making computational cloning design accessible through natural language interaction with AI assistants.*

## About MCP (Model Context Protocol)

MCP is a protocol that bridges the gap between AI systems and specialized domain knowledge. For molecular biology, it enables:

- **Structured DNA Operations**: Direct access to PyDNA's simulation capabilities
- **Natural Language Cloning**: Describe experiments in plain English, get computational results
- **Type Safety**: Strong typing and validation for molecular biology data
- **AI Integration**: Seamless integration with AI assistants for experimental design

## Available Tools

The PyDNA MCP server provides comprehensive tools organized by functionality:

### 🧬 Sequence Operations
- **`pydna_create_sequence`** - Create DNA sequences from strings with optional circularity
- **`pydna_sequence_info`** - Get detailed sequence analysis (GC content, Tm, MW, etc.)
- **`pydna_reverse_complement`** - Generate reverse complement of DNA sequences

### 🔬 PCR & Primer Design  
- **`pydna_design_primers`** - Automated primer design with customizable melting temperatures
- **`pydna_pcr_amplify`** - Simulate PCR amplification with primers and template
- **`pydna_anneal_primers`** - Analyze primer annealing to templates

### ✂️ Restriction Analysis
- **`pydna_restriction_analysis`** - Map restriction enzyme cutting sites
- **`pydna_digest_sequence`** - Simulate restriction digestion and fragment analysis

### 🔧 DNA Assembly
- **`pydna_assembly`** - Homologous recombination-based assembly
- **`pydna_gibson_assembly`** - Gibson assembly design and simulation

### 📁 File Operations
- **`pydna_read_sequence_file`** - Read FASTA/GenBank/other sequence files
- **`pydna_write_sequence_file`** - Write sequences to various file formats
- **`pydna_download_genbank`** - Download sequences directly from GenBank

### ⚡ Analysis Tools
- **`pydna_gel_electrophoresis`** - Simulate gel electrophoresis runs

## Available Resources

- **`resource://pydna_help`** - Comprehensive help with workflows and examples
- **`resource://pydna_examples`** - Common molecular cloning workflow examples

## Quick Start

### Installing uv

```bash
# Download and install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Verify installation
uv --version
uvx --version
```

uvx is a very nice tool that can run a python package installing it if needed.

### Running with uvx

You can run the pydna-mcp server directly using uvx without cloning the repository:

```bash
# Run the server in STDIO mode (default)
uvx pydna-mcp
```

<details>
<summary>Other uvx modes (STDIO, HTTP, SSE)</summary>

#### STDIO Mode (for MCP clients that require stdio)

```bash
# Or explicitly specify stdio mode  
uvx pydna-mcp stdio
```

#### HTTP Mode (Web Server)
```bash
# Run the server in streamable HTTP mode on default (3001) port
uvx pydna-mcp server

# Run on a specific port
uvx pydna-mcp server --port 8000
```

#### SSE Mode (Server-Sent Events)
```bash
# Run the server in SSE mode
uvx pydna-mcp sse
```

</details>

In cases when there are problems with uvx often they can be caused by cleaning uv cache:
```
uv cache clean
```

The HTTP mode will start a web server that you can access at `http://localhost:3001/mcp` (with documentation at `http://localhost:3001/docs`). The STDIO mode is designed for MCP clients that communicate via standard input/output, while SSE mode uses Server-Sent Events for real-time communication.

## Configuring your AI Client (Anthropic Claude Desktop, Cursor, Windsurf, etc.)

### Quick Configuration Example

Here's what you can copy directly into your Claude Desktop or Cursor MCP configuration:

```json
{
  "mcpServers": {
    "pydna-mcp": {
      "command": "uvx",
      "args": ["pydna-mcp"],
      "env": {
        "MCP_PORT": "3001",
        "MCP_HOST": "0.0.0.0", 
        "MCP_TRANSPORT": "stdio"
      }
    }
  }
}
```

### Alternative: Using Preconfigured Files

We also provide preconfigured JSON files for different use cases (coming soon):

- **For STDIO mode (recommended):** Use `mcp-config-stdio.json`
- **For HTTP mode:** Use `mcp-config.json`
- **For local development:** Use `mcp-config-stdio-debug.json`

### Inspecting PyDNA MCP server

<details>
<summary>Using MCP Inspector to explore server capabilities</summary>

If you want to inspect the methods provided by the MCP server, use npx (you may need to install nodejs and npm):

For STDIO mode with uvx:
```bash
npx @modelcontextprotocol/inspector uvx pydna-mcp
```

You can also run the inspector manually and configure it through the interface:
```bash
npx @modelcontextprotocol/inspector
```

After that you can explore the tools and resources with MCP Inspector at http://127.0.0.1:6274 (note, if you run inspector several times it can change port)

</details>

### Integration with AI Systems

Simply point your AI client (like Cursor, Windsurf, ClaudeDesktop, VS Code with Copilot, or [others](https://github.com/punkpeye/awesome-mcp-clients)) to use the appropriate configuration from above.

## Repository setup

```bash
# Clone the repository
git clone https://github.com/antonkulaga/pydna-mcp.git
cd pydna-mcp
uv sync
```

### Running the MCP Server

If you already cloned the repo you can run the server with uv:

```bash
# Start the MCP server locally (HTTP mode)
uv run server

# Or start in STDIO mode
uv run stdio

# Or start in SSE mode  
uv run sse
```

## Common Molecular Biology Workflows

<details>
<summary>Detailed workflow examples</summary>

### 1. Basic PCR Design & Simulation

```python
# Natural language workflow:
# "Design primers for my gene and simulate PCR"

# 1. Create your template sequence
template = "ATGGCTAGCATGACTGGTGGACAGCAAATGGGTCG..."

# 2. Design primers automatically
primers = design_primers(template, target_tm=60.0)

# 3. Simulate PCR amplification  
amplicon = pcr_amplify(primers.forward, primers.reverse, template)

# 4. Analyze the product
gel_result = gel_electrophoresis([amplicon.sequence])
```

### 2. Restriction Cloning Workflow

```python
# Natural language workflow:
# "I want to clone my PCR product into pUC19 using EcoRI and BamHI"

# 1. Analyze restriction sites in your insert
insert_analysis = restriction_analysis(pcr_product, ["EcoRI", "BamHI"])

# 2. Analyze sites in vector
vector_analysis = restriction_analysis(puc19_sequence, ["EcoRI", "BamHI"])  

# 3. Digest both sequences
digested_insert = digest_sequence(pcr_product, ["EcoRI", "BamHI"])
digested_vector = digest_sequence(puc19_sequence, ["EcoRI", "BamHI"])

# 4. Check fragment sizes on gel
gel_check = gel_electrophoresis([digested_insert.fragments, digested_vector.fragments])
```

### 3. Gibson Assembly Design

```python
# Natural language workflow:  
# "Design Gibson assembly to join three DNA fragments"

# 1. Prepare your fragments with overlaps
fragments = [fragment1, fragment2, fragment3]

# 2. Design Gibson assembly
assembly_result = gibson_assembly(fragments, overlap=20)

# 3. Verify the final construct
final_info = sequence_info(assembly_result.sequence)
```

### 4. Primer Analysis & Optimization

```python
# Natural language workflow:
# "Check if my existing primers will work with this template"

# 1. Test primer annealing
annealing_result = anneal_primers([forward_primer, reverse_primer], template)

# 2. If primers need optimization, redesign them
if annealing_result.efficiency < 0.8:
    better_primers = design_primers(template, target_tm=58.0)
```

</details>

## Example Use Cases

<details>
<summary>Research applications you can accomplish with PyDNA MCP</summary>

### Gene Cloning & Expression
- Design primers for gene amplification from genomic DNA
- Plan restriction-based cloning into expression vectors
- Optimize primer sequences for difficult templates
- Design Gibson assembly for multi-fragment constructs

### Plasmid Construction
- Analyze existing plasmid maps and restriction sites
- Plan subcloning strategies between different vectors  
- Design seamless assembly workflows
- Verify construct integrity through simulation

### PCR Optimization
- Test different primer combinations computationally
- Analyze primer-dimer formation potential
- Optimize annealing temperatures for multiplex PCR
- Design primers with specific restriction sites

### Molecular Biology Education
- Simulate lab experiments before bench work
- Understand restriction mapping concepts
- Learn PCR design principles through practice
- Explore DNA assembly mechanisms

### Synthetic Biology  
- Design modular DNA parts with standardized overlaps
- Plan complex multi-fragment assemblies
- Optimize construct sequences for expression
- Simulate synthetic pathway construction

</details>

## Safety Features

- **Computational Only**: All operations are simulation-based, no wet lab risks
- **Input Validation**: Comprehensive validation of DNA sequences and parameters  
- **Error Handling**: Clear error messages for invalid sequences or operations
- **Type Safety**: Strongly typed interfaces prevent common errors

## Testing & Development

The MCP server includes comprehensive tests for all molecular biology operations:

### Running Tests

```bash
# Run all tests  
uv run pytest -vvv -s

# Run specific test categories
uv run pytest tests/test_pydna_mcp_server.py -v
```

### Test Coverage

Tests cover:
- Sequence creation and manipulation
- PCR simulation accuracy
- Primer design algorithms  
- Restriction analysis correctness
- Assembly simulation results
- File format compatibility
- Error handling scenarios

## Contributing

We welcome contributions from the molecular biology and bioinformatics community! 🧬

Whether you're a researcher, developer, or student interested in computational molecular biology, there are many ways to contribute:

### Ways to Contribute

- **🐛 Bug Reports**: Found issues with simulations? Report them!
- **💡 Feature Requests**: Need additional PyDNA functionality? Let us know!
- **📝 Documentation**: Improve examples, tutorials, or workflow guides
- **🧪 Testing**: Add test cases for molecular biology edge cases
- **🔬 Validation**: Compare simulation results with experimental data
- **⚡ Performance**: Optimize calculations for large sequences
- **🌐 Integration**: Create examples for new MCP clients
- **🎥 Tutorials**: Create videos showing molecular cloning workflows
- **📚 Educational Content**: Develop teaching materials using PyDNA MCP

**Tutorials and educational content are especially valuable!** Computational molecular biology tools can significantly accelerate research and education when properly demonstrated.

### Getting Started

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-cloning-tool`)
3. Make your changes and add tests
4. Run the test suite (`uv run pytest`)  
5. Commit your changes (`git commit -m 'Add amazing cloning tool'`)
6. Push to your branch (`git push origin feature/amazing-cloning-tool`)
7. Open a Pull Request

### Development Guidelines

- Follow existing code style (we use `black` for formatting)
- Add tests for new molecular biology functionality
- Update documentation for new features
- Validate biological accuracy of simulations
- Write clear commit messages

## Known Issues

### PyDNA Library Limitations
This MCP server exposes PyDNA's capabilities, which means it inherits any limitations of the underlying library. Some advanced cloning techniques or specialized molecular biology applications may require additional tool development.

### Sequence Size Limitations  
Very large sequences (>100kb) may experience performance issues. Optimization for genomic-scale sequences is ongoing.

## License

This project is licensed under the MIT License.

## Acknowledgments

- **[PyDNA](https://github.com/BjornFJohansson/pydna)** for the comprehensive DNA simulation library
  - Pereira, F., Azevedo, F., Carvalho, Â., Ribeiro, G. F., Budde, M. W., & Johansson, B. (2015). Pydna: a simulation and documentation tool for DNA assembly strategies using python. BMC bioinformatics, 16(1), 142.
- **[BioPython](https://biopython.org/)** for foundational biological sequence analysis
- **[Model Context Protocol](https://modelcontextprotocol.io/)** for the protocol specification  
- **[FastMCP](https://github.com/jlowin/fastmcp)** for the MCP server framework

### Other MCP Servers in the Holy Bio MCP Framework

This server is part of the complete **[Holy Bio MCP](https://github.com/longevity-genie/holy-bio-mcp)** framework, which includes:

- **[gget-mcp](https://github.com/longevity-genie/gget-mcp)** - Powerful bioinformatics toolkit for genomics queries and analysis
- **[opengenes-mcp](https://github.com/longevity-genie/opengenes-mcp)** - Database of genes involved in aging and longevity
- **[synergy-age-mcp](https://github.com/longevity-genie/synergy-age-mcp)** - Database of synergistic and antagonistic genetic interactions in longevity
- **[biothings-mcp](https://github.com/longevity-genie/biothings-mcp)** - Access to BioThings.io APIs for comprehensive gene, variant, chemical, and taxonomic data
- **[pharmacology-mcp](https://github.com/antonkulaga/pharmacology-mcp)** - Access to the Guide to PHARMACOLOGY database

The framework provides unified configuration files that enable all servers at once, making it easy to access 50+ specialized bioinformatics functions through a single setup. This award-winning project continues to evolve as a comprehensive platform for AI-driven biological research.
