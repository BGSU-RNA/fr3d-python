"""
Pydantic models for R3DCID (Circular Interaction Diagrams) form inputs.

This module defines abstract base classes and Pydantic input field models
for validating and parsing R3DCID form submissions and URL parameters.
"""

from abc import ABC, abstractmethod
from enum import Enum
from typing import Annotated, Optional

from pydantic import BaseModel, Field, field_validator, model_validator


# ============================================================================
# Enums for form options
# ============================================================================

class OutputFormat(str, Enum):
    """Output format options for the circular diagram."""
    HTML = "html"
    SVG = "svg"
    PDF = "pdf"


class ColorScheme(str, Enum):
    """Color scheme options for the diagram."""
    DEFAULT = "default"
    WONG = "wong"  # Colorblind-safe palette
    GRAYSCALE = "grayscale"


class InteractionVisibility(str, Enum):
    """Visibility options for interaction arcs."""
    SHOW = "show"
    DIM = "dim"
    HIDE = "hide"


class TextOption(str, Enum):
    """Options for text displayed outside the circle."""
    HELIX = "helix"
    BASEPAIR = "basepair"
    STACKING = "stacking"
    BPH = "bph"  # Base-phosphate
    BR = "br"   # Base-ribose
    SR = "sr"   # Sugar-ribose
    SO = "so"   # Oxygen stacking
    NEAR = "near"
    ALL = "all"


class HeaderField(str, Enum):
    """Header fields to display in output."""
    TITLE = "title"
    METHOD = "method"
    RELEASE_DATE = "release_date"
    SOURCE = "source"
    RESOLUTION = "resolution"


# ============================================================================
# Interaction Arc Settings Model
# ============================================================================

class InteractionArcSettings(BaseModel):
    """Settings for interaction arc visibility.
    
    Each interaction type can be set to show, dim, or hide.
    All default to 'show'.
    """
    nested_wc: InteractionVisibility = Field(
        default=InteractionVisibility.SHOW,
        description="Nested Watson-Crick basepairs visibility"
    )
    lr_wc: InteractionVisibility = Field(
        default=InteractionVisibility.SHOW,
        description="Long-range Watson-Crick basepairs visibility"
    )
    nested_non_wc: InteractionVisibility = Field(
        default=InteractionVisibility.SHOW,
        description="Nested non-Watson-Crick basepairs visibility"
    )
    lr_non_wc: InteractionVisibility = Field(
        default=InteractionVisibility.SHOW,
        description="Long-range non-Watson-Crick basepairs visibility"
    )
    stacking: InteractionVisibility = Field(
        default=InteractionVisibility.SHOW,
        description="Base stacking interactions visibility"
    )
    bph: InteractionVisibility = Field(
        default=InteractionVisibility.SHOW,
        description="Base-phosphate interactions visibility"
    )
    br: InteractionVisibility = Field(
        default=InteractionVisibility.SHOW,
        description="Base-ribose interactions visibility"
    )
    sr: InteractionVisibility = Field(
        default=InteractionVisibility.SHOW,
        description="Sugar-ribose interactions visibility"
    )
    so: InteractionVisibility = Field(
        default=InteractionVisibility.SHOW,
        description="Oxygen stacking interactions visibility"
    )
    near: InteractionVisibility = Field(
        default=InteractionVisibility.SHOW,
        description="Near interactions visibility"
    )

    def get_show_list(self) -> list[str]:
        """Return list of interaction types set to show."""
        result = []
        for field_name, value in self:
            if value == InteractionVisibility.SHOW:
                result.append(field_name.replace("_", "-"))
        return result

    def get_dim_list(self) -> list[str]:
        """Return list of interaction types set to dim."""
        result = []
        for field_name, value in self:
            if value == InteractionVisibility.DIM:
                result.append(field_name.replace("_", "-"))
        return result

    def get_hide_list(self) -> list[str]:
        """Return list of interaction types set to hide."""
        result = []
        for field_name, value in self:
            if value == InteractionVisibility.HIDE:
                result.append(field_name.replace("_", "-"))
        return result


# ============================================================================
# Main R3DCID Input Model
# ============================================================================

class R3DCIDInput(BaseModel):
    """
    Pydantic model for R3DCID form input validation.
    
    This model validates all input parameters for generating circular
    interaction diagrams for RNA and DNA 3D structures.
    
    Example usage:
        input_data = R3DCIDInput(
            chains="7EZ2",
            output_format=OutputFormat.HTML,
            color_scheme=ColorScheme.DEFAULT
        )
    """
    
    # Required fields
    chains: Annotated[
        str,
        Field(
            min_length=4,
            max_length=500,
            description=(
                "PDB id, models, chains specification. "
                "Format: 'pdbid|model1|chain1+pdbid|model2|chain2+...'"
            ),
            examples=["7EZ2", "7K00", "8GLP|1|L8+8GLP|1|L5"],
        )
    ]
    
    # Optional structure specification fields
    assemblies: Annotated[
        Optional[str],
        Field(
            default=None,
            max_length=100,
            description="Comma-separated assembly numbers (e.g., '4,2')",
            examples=["4,2", "1"],
        )
    ] = None
    
    symmetries: Annotated[
        Optional[str],
        Field(
            default=None,
            max_length=100,
            description="Comma-separated symmetry operators",
        )
    ] = None
    
    # Output format options
    output_format: OutputFormat = Field(
        default=OutputFormat.HTML,
        description="Output format for the diagram"
    )
    
    description: Annotated[
        Optional[str],
        Field(
            default=None,
            max_length=300,
            description=(
                "Custom description for the header. "
                "Use \\n for line breaks."
            ),
        )
    ] = None
    
    # Coloring options
    color_scheme: ColorScheme = Field(
        default=ColorScheme.DEFAULT,
        alias="coloring",
        description="Color scheme for the diagram"
    )
    
    # Interaction arc visibility
    interaction_arcs: InteractionArcSettings = Field(
        default_factory=InteractionArcSettings,
        description="Visibility settings for interaction arcs"
    )
    
    # Text options (checkboxes)
    text_options: list[TextOption] = Field(
        default=[TextOption.BASEPAIR],
        description="Types of text to display outside the circle"
    )
    
    # Helix size override
    helix_size: Annotated[
        Optional[int],
        Field(
            default=None,
            ge=0,
            le=100,
            alias="hs",
            description=(
                "Override font size for helix numbers. "
                "0 means no helix numbers."
            ),
        )
    ] = None
    
    # Display nucleotides without 3D coordinates
    show_no_3d_coords: Annotated[
        bool,
        Field(
            default=True,
            alias="n3d",
            description="Whether to display nucleotides without 3D coordinates"
        )
    ] = True
    
    # Header fields to display
    header_fields: list[HeaderField] = Field(
        default=[
            HeaderField.TITLE,
            HeaderField.METHOD,
            HeaderField.RELEASE_DATE,
            HeaderField.SOURCE,
            HeaderField.RESOLUTION,
        ],
        description="Fields to display in the output header"
    )

    model_config = {
        "populate_by_name": True,
        "use_enum_values": True,
        "str_strip_whitespace": True,
    }

    @field_validator("chains")
    @classmethod
    def validate_chains(cls, v: str) -> str:
        """Validate and normalize chains input."""
        # Remove extra whitespace and normalize separators
        v = v.strip().replace(" ", "+")
        
        # Basic validation - should contain at least a PDB ID (4 chars)
        if not v or len(v) < 4:
            raise ValueError("Chains must contain at least a valid PDB ID (4 characters)")
        
        return v

    @field_validator("assemblies", "symmetries")
    @classmethod
    def validate_comma_separated(cls, v: Optional[str]) -> Optional[str]:
        """Validate and normalize comma-separated values."""
        if v is None or v.strip() == "":
            return None
        # Normalize separators
        return v.strip().replace(";", ",")

    @model_validator(mode="after")
    def validate_text_options(self) -> "R3DCIDInput":
        """If 'all' is in text_options, expand to all options."""
        if TextOption.ALL in self.text_options:
            self.text_options = [
                TextOption.HELIX,
                TextOption.BASEPAIR,
                TextOption.STACKING,
                TextOption.BPH,
                TextOption.BR,
                TextOption.SR,
                TextOption.SO,
                TextOption.NEAR,
            ]
        return self

    def to_url_params(self) -> dict[str, str]:
        """Convert model to URL parameters dict."""
        params: dict[str, str] = {"chains": self.chains}
        
        if self.assemblies:
            params["assemblies"] = self.assemblies
        
        if self.symmetries:
            params["symmetries"] = self.symmetries
        
        # Handle output_format (may be enum or string due to use_enum_values)
        fmt = self.output_format.value if hasattr(self.output_format, 'value') else self.output_format
        if fmt != "pdf":
            params["format"] = fmt
        
        # Handle color_scheme (may be enum or string due to use_enum_values)
        color = self.color_scheme.value if hasattr(self.color_scheme, 'value') else self.color_scheme
        if color != "default":
            params["coloring"] = color
        
        if self.description:
            params["description"] = self.description
        
        # Get hide and dim lists from interaction arcs
        hide_list = self.interaction_arcs.get_hide_list()
        dim_list = self.interaction_arcs.get_dim_list()
        
        if hide_list:
            params["hide"] = ",".join(hide_list)
        
        if dim_list:
            params["dim"] = ",".join(dim_list)
        
        # Text options (may be enum or string due to use_enum_values)
        text_values = [
            opt.value if hasattr(opt, 'value') else opt 
            for opt in self.text_options
        ]
        if text_values != ["basepair"]:
            params["text"] = ",".join(text_values)
        
        if self.helix_size is not None:
            params["hs"] = str(self.helix_size)
        
        if not self.show_no_3d_coords:
            params["n3d"] = "false"
        
        # Header fields (may be enum or string due to use_enum_values)
        header_values = [
            f.value if hasattr(f, 'value') else f 
            for f in self.header_fields
        ]
        default_headers = ["title", "method", "release_date", "source", "resolution"]
        if sorted(header_values) != sorted(default_headers):
            if header_values:
                params["header"] = ",".join(header_values)
            else:
                params["header"] = "none"
        
        return params


# ============================================================================
# Request model for parsing URL query parameters
# ============================================================================

class R3DCIDQueryParams(BaseModel):
    """
    Model for parsing R3DCID query parameters from URL.
    
    This model handles the flat query parameter format used in URLs,
    converting them to the structured R3DCIDInput format.
    """
    
    chains: str = Field(
        ...,
        min_length=4,
        description="PDB id and chain specification"
    )
    assemblies: Optional[str] = Field(default=None)
    symmetries: Optional[str] = Field(default=None)
    format: Optional[str] = Field(default="html")
    coloring: Optional[str] = Field(default="default")
    description: Optional[str] = Field(default=None, max_length=300)
    show: Optional[str] = Field(default=None)
    dim: Optional[str] = Field(default=None)
    hide: Optional[str] = Field(default=None)
    text: Optional[str] = Field(default=None)
    hs: Optional[str] = Field(default=None)
    n3d: Optional[str] = Field(default="true")
    header: Optional[str] = Field(default=None)
    input_form: Optional[str] = Field(default=None)

    model_config = {
        "str_strip_whitespace": True,
    }

    def to_r3dcid_input(self) -> R3DCIDInput:
        """Convert query parameters to R3DCIDInput model."""
        # Parse format
        output_format = OutputFormat.HTML
        if self.format:
            try:
                output_format = OutputFormat(self.format.lower())
            except ValueError:
                output_format = OutputFormat.PDF
        
        # Parse coloring
        color_scheme = ColorScheme.DEFAULT
        if self.coloring:
            try:
                color_scheme = ColorScheme(self.coloring.lower())
            except ValueError:
                color_scheme = ColorScheme.DEFAULT
        
        # Parse interaction arc visibility
        arc_settings = InteractionArcSettings()
        arc_field_map = {
            "nested-wc": "nested_wc",
            "lr-wc": "lr_wc",
            "nested-non-wc": "nested_non_wc",
            "lr-non-wc": "lr_non_wc",
            "stacking": "stacking",
            "bph": "bph",
            "br": "br",
            "sr": "sr",
            "so": "so",
            "near": "near",
        }
        
        if self.hide:
            for item in self.hide.split(","):
                item = item.strip().lower()
                if item in arc_field_map:
                    setattr(arc_settings, arc_field_map[item], InteractionVisibility.HIDE)
        
        if self.dim:
            for item in self.dim.split(","):
                item = item.strip().lower()
                if item in arc_field_map:
                    setattr(arc_settings, arc_field_map[item], InteractionVisibility.DIM)
        
        if self.show:
            # If show is specified, first set all to hide, then show specified ones
            for field_name in arc_field_map.values():
                setattr(arc_settings, field_name, InteractionVisibility.HIDE)
            for item in self.show.split(","):
                item = item.strip().lower()
                if item in arc_field_map:
                    setattr(arc_settings, arc_field_map[item], InteractionVisibility.SHOW)
        
        # Parse text options
        text_options: list[TextOption] = []
        if self.text:
            for item in self.text.split(","):
                item = item.strip().lower()
                try:
                    text_options.append(TextOption(item))
                except ValueError:
                    pass
        if not text_options:
            text_options = [TextOption.BASEPAIR]
        
        # Parse helix size
        helix_size: Optional[int] = None
        if self.hs:
            try:
                helix_size = int(self.hs)
            except ValueError:
                pass
        
        # Parse n3d
        show_no_3d = True
        if self.n3d and self.n3d.lower() == "false":
            show_no_3d = False
        
        # Parse header fields
        header_fields: list[HeaderField] = []
        if self.header:
            if self.header.lower() == "none":
                header_fields = []
            else:
                header_map = {
                    "title": HeaderField.TITLE,
                    "method": HeaderField.METHOD,
                    "release_date": HeaderField.RELEASE_DATE,
                    "release": HeaderField.RELEASE_DATE,
                    "source": HeaderField.SOURCE,
                    "resolution": HeaderField.RESOLUTION,
                }
                for item in self.header.split(","):
                    item = item.strip().lower()
                    if item in header_map:
                        header_fields.append(header_map[item])
        else:
            header_fields = [
                HeaderField.TITLE,
                HeaderField.METHOD,
                HeaderField.RELEASE_DATE,
                HeaderField.SOURCE,
                HeaderField.RESOLUTION,
            ]
        
        return R3DCIDInput(
            chains=self.chains,
            assemblies=self.assemblies,
            symmetries=self.symmetries,
            output_format=output_format,
            color_scheme=color_scheme,
            description=self.description,
            interaction_arcs=arc_settings,
            text_options=text_options,
            helix_size=helix_size,
            show_no_3d_coords=show_no_3d,
            header_fields=header_fields,
        )


# ============================================================================
# Data Models for Structure Information
# ============================================================================

class ChainInfo(BaseModel):
    """Information about a chain in a PDB structure."""
    pdb_id: str = Field(description="PDB identifier")
    chain_name: str = Field(description="Chain identifier (e.g., 'A', 'L5')")
    model: str = Field(default="1", description="Model number")
    assembly: str = Field(default="1", description="Assembly identifier")
    symmetry: Optional[str] = Field(default=None, description="Symmetry operator")
    entity_type: Optional[str] = Field(
        default=None,
        description="Entity type (rna, dna, hybrid)"
    )
    rfam: Optional[str] = Field(default=None, description="Rfam family identifier")
    chain_length: int = Field(default=0, description="Number of nucleotides")
    chain_priority: int = Field(default=0, description="Priority for ordering")
    symmetry_priority: int = Field(default=0, description="Symmetry priority")
    assembly_priority: int = Field(default=0, description="Assembly priority")
    model_chain_priority: int = Field(default=0, description="Model-chain priority")
    final_cww_group: float = Field(
        default=0.0,
        description="Final cWW grouping for ordering around diagram"
    )
    
    model_config = {"extra": "allow"}


class AssemblyInfo(BaseModel):
    """Information about assemblies in a structure."""
    ok_assemblies: list[str] = Field(
        default_factory=list,
        description="List of valid assembly identifiers"
    )
    valid_assembly_pairs: list[tuple[str, str]] = Field(
        default_factory=list,
        description="List of valid (assembly1, assembly2) pairs for interactions"
    )
    message: str = Field(default="", description="Status or error message")


class HeaderInfo(BaseModel):
    """PDB structure header information."""
    title: Optional[str] = Field(default=None, description="Structure title")
    method: Optional[str] = Field(default=None, description="Experimental method")
    release_date: Optional[str] = Field(default=None, description="Release date")
    source: Optional[str] = Field(default=None, description="Biological source organism")
    resolution: Optional[str] = Field(default=None, description="Resolution in Angstroms")


class DiagramOutput(BaseModel):
    """Output from diagram generation."""
    filename: str = Field(description="Base filename (without extension)")
    ps_commands: str = Field(default="", description="PostScript commands")
    svg_commands: str = Field(default="", description="SVG commands")
    page_width: float = Field(default=612, description="Page width in points")
    page_height: float = Field(default=792, description="Page height in points")
    output_path: Optional[str] = Field(default=None, description="Full output path")
    message: str = Field(default="", description="Status or error message")


class InteractionTriple(BaseModel):
    """An interaction between two nucleotides with crossing count."""
    unit_id_1: str = Field(description="First nucleotide unit ID")
    unit_id_2: str = Field(description="Second nucleotide unit ID")
    crossing: int = Field(
        default=0,
        description="Number of nested cWW basepairs crossed"
    )


class SequenceMapping(BaseModel):
    """Mapping from sequence position to unit ID."""
    sequence_id: str = Field(description="Sequence position identifier")
    unit_id: str = Field(description="Unit ID (or 'NULL' if no 3D coordinates)")
    
    @property
    def has_3d_coords(self) -> bool:
        """Check if this position has 3D coordinates."""
        return self.unit_id != "NULL"


# ============================================================================
# Abstract Base Classes for Users
# ============================================================================

class AbstractDiagramGenerator(ABC):
    """Abstract base class for circular diagram generators."""
    
    @abstractmethod
    def generate(self, input_data: R3DCIDInput) -> DiagramOutput:
        """
        Generate a circular interaction diagram.
        
        Args:
            input_data: Validated R3DCID input parameters.
            
        Returns:
            DiagramOutput with paths and generated content.
        """
        pass
    
    @abstractmethod
    def get_filename(self, input_data: R3DCIDInput) -> str:
        """
        Generate the output filename based on input parameters.
        
        Args:
            input_data: Validated R3DCID input parameters.
            
        Returns:
            Base filename (without extension).
        """
        pass
    
    @abstractmethod
    def draw_circular_diagram(
        self,
        chain_info: list[ChainInfo],
        assemblies: AssemblyInfo,
        filename: str,
        interaction_to_triple_list: dict[str, list[InteractionTriple]],
        params: dict,
    ) -> tuple[str, str, float, float]:
        """
        Construct PostScript and SVG strings for the circular diagram.
        
        Args:
            chain_info: List of chain information dictionaries.
            assemblies: Assembly information.
            filename: Base filename.
            interaction_to_triple_list: Interactions by type.
            params: Display parameters.
            
        Returns:
            Tuple of (PS commands, SVG commands, page_width, page_height).
        """
        pass
    
    @abstractmethod
    def draw_arcs(
        self,
        pairs_and_crossing: list[InteractionTriple],
        unit_id_to_angle: dict[str, float],
        arc_group: str,
        params: dict,
    ) -> tuple[str, str]:
        """
        Draw interaction arcs for a specific interaction type.
        
        Args:
            pairs_and_crossing: List of interaction triples.
            unit_id_to_angle: Mapping from unit ID to angle on circle.
            arc_group: The arc group name (e.g., 'nested-wc', 'stacking').
            params: Display parameters.
            
        Returns:
            Tuple of (PS commands, SVG commands).
        """
        pass


class AbstractStructureReader(ABC):
    """Abstract base class for reading structure data."""
    
    @abstractmethod
    def read_structure(self, pdb_id: str) -> dict:
        """
        Read structure data for a PDB ID.
        
        Args:
            pdb_id: The PDB identifier.
            
        Returns:
            Dictionary containing structure data.
        """
        pass
    
    @abstractmethod
    def get_chains(self, pdb_id: str) -> list[str]:
        """
        Get list of available chains for a structure.
        
        Args:
            pdb_id: The PDB identifier.
            
        Returns:
            List of chain identifiers.
        """
        pass
    
    @abstractmethod
    def get_chain_info(self, pdb_id: str) -> list[ChainInfo]:
        """
        Get detailed chain information from the RNA 3D Hub API.
        
        Args:
            pdb_id: The PDB identifier.
            
        Returns:
            List of ChainInfo for all chains in the structure.
        """
        pass
    
    @abstractmethod
    def get_header_info(self, pdb_id: str) -> HeaderInfo:
        """
        Get PDB header information (title, method, resolution, etc.).
        
        Args:
            pdb_id: The PDB identifier.
            
        Returns:
            HeaderInfo with structure metadata.
        """
        pass
    
    @abstractmethod
    def get_sequence_mapping(
        self,
        pdb_id: str,
        chain: str,
        model: str = "1"
    ) -> list[SequenceMapping]:
        """
        Get sequence position to unit ID mapping for a chain.
        
        Args:
            pdb_id: The PDB identifier.
            chain: Chain identifier.
            model: Model number (default "1").
            
        Returns:
            List of SequenceMapping for the chain.
        """
        pass
    
    @abstractmethod
    def get_nucleotide_annotations(
        self,
        pdb_id: str,
        chain: str,
        model: str = "1"
    ) -> dict[str, str]:
        """
        Get nucleotide annotations (e.g., secondary structure) for a chain.
        
        Args:
            pdb_id: The PDB identifier.
            chain: Chain identifier.
            model: Model number (default "1").
            
        Returns:
            Dictionary mapping unit_id to annotation string.
        """
        pass


class AbstractInteractionReader(ABC):
    """Abstract base class for reading interaction annotations."""
    
    @abstractmethod
    def read_interactions(self, pdb_id: str) -> dict[str, list[InteractionTriple]]:
        """
        Read FR3D interaction annotations for a structure.
        
        Args:
            pdb_id: The PDB identifier.
            
        Returns:
            Dictionary mapping interaction types to lists of InteractionTriple.
        """
        pass
    
    @abstractmethod
    def download_interactions(
        self,
        pdb_id: str,
        data_directory: str = ""
    ) -> tuple[dict[str, list], Optional[str]]:
        """
        Download and cache FR3D interaction annotations.
        
        Args:
            pdb_id: The PDB identifier.
            data_directory: Directory for caching files.
            
        Returns:
            Tuple of (interaction dict, error message or None).
        """
        pass
    
    @abstractmethod
    def get_interaction_types(self) -> list[str]:
        """
        Get list of all supported interaction types.
        
        Returns:
            List of interaction type strings (e.g., 'cWW', 's35', '5BPh').
        """
        pass


class AbstractChainOrderer(ABC):
    """Abstract base class for ordering chains around the diagram."""
    
    @abstractmethod
    def order_chains(
        self,
        pdb_id: str,
        requested_assemblies: list[str],
        requested_models: list[str],
        requested_model_chain: list[tuple[str, str]],
        requested_chains: list[str],
        requested_symmetries: list[str],
        interaction_to_triple_list: dict[str, list],
    ) -> tuple[list[ChainInfo], AssemblyInfo, Optional[str]]:
        """
        Determine the order of chains around the circular diagram.
        
        Args:
            pdb_id: The PDB identifier.
            requested_assemblies: List of requested assembly IDs.
            requested_models: List of requested model numbers.
            requested_model_chain: List of (model, chain) tuples.
            requested_chains: List of requested chain IDs.
            requested_symmetries: List of requested symmetry operators.
            interaction_to_triple_list: Interactions by type.
            
        Returns:
            Tuple of (sorted chain info list, assembly info, error message).
        """
        pass
    
    @abstractmethod
    def group_chains_by_cww(
        self,
        chain_info: list[ChainInfo],
        interaction_to_triple_list: dict[str, list],
    ) -> list[ChainInfo]:
        """
        Group chains that have cWW basepairs between them.
        
        Args:
            chain_info: List of chain information.
            interaction_to_triple_list: Interactions by type.
            
        Returns:
            Chain info list with updated cWW group assignments.
        """
        pass


class AbstractInputProcessor(ABC):
    """Abstract base class for processing user input."""
    
    @abstractmethod
    def process_input_chains(
        self,
        input_text: str,
        params: Optional[dict] = None,
    ) -> tuple[str, str, list[str], list[str], list[tuple[str, str]], list[str], list[str]]:
        """
        Parse input text specifying PDB id, models, chains.
        
        Args:
            input_text: User input string like '8GLP|1|L5+8GLP|1|L8'.
            params: Optional parameters dict with 'assembly', 'symmetry'.
            
        Returns:
            Tuple of:
                - pdb_id
                - filename
                - requested_assemblies
                - requested_models  
                - requested_model_chain (list of tuples)
                - requested_chains
                - requested_symmetries
        """
        pass
    
    @abstractmethod
    def set_parameters_from_input(
        self,
        params: dict,
        filename: str,
        pdb_id: str,
    ) -> tuple[dict, str]:
        """
        Process display parameters and update filename accordingly.
        
        Args:
            params: Input parameters dictionary.
            filename: Current filename.
            pdb_id: PDB identifier.
            
        Returns:
            Tuple of (updated params, updated filename).
        """
        pass



# ============================================================================
# URL Builder Functions
# ============================================================================

DEFAULT_BASE_URL = "https://rna.bgsu.edu/fr3d/r3dcid"


def build_r3dcid_url(
    chains: str,
    *,
    assemblies: Optional[str] = None,
    symmetries: Optional[str] = None,
    output_format: OutputFormat | str = OutputFormat.HTML,
    coloring: ColorScheme | str = ColorScheme.DEFAULT,
    description: Optional[str] = None,
    show: Optional[list[str] | str] = None,
    dim: Optional[list[str] | str] = None,
    hide: Optional[list[str] | str] = None,
    text: Optional[list[str] | str] = None,
    helix_size: Optional[int] = None,
    n3d: bool = True,
    header: Optional[list[str] | str] = None,
    input_form: bool = False,
    base_url: str = DEFAULT_BASE_URL,
) -> str:
    """
    Build an R3DCID URL with the specified parameters.
    
    Args:
        chains: PDB id, models, chains specification.
                Format: 'pdbid|model|chain+pdbid|model|chain+...'
                Examples: '7EZ2', '8GLP|1|L8+8GLP|1|L5', '7K00'
        assemblies: Comma-separated assembly numbers (e.g., '4,2')
        symmetries: Comma-separated symmetry operators
        output_format: Output format ('html', 'svg', 'pdf') or OutputFormat enum
        coloring: Color scheme ('default', 'wong', 'grayscale') or ColorScheme enum
        description: Custom description for the header (use \\n for line breaks)
        show: Interaction types to show (others hidden).
              Options: nested-wc, lr-wc, nested-non-wc, lr-non-wc, stacking, bph, br, sr, so, near
        dim: Interaction types to dim.
        hide: Interaction types to hide.
        text: Text options to display outside circle.
              Options: helix, basepair, stacking, bph, br, sr, so, near, all
        helix_size: Override font size for helix numbers (0 = no helix numbers)
        n3d: Whether to display nucleotides without 3D coordinates (default True)
        header: Header fields to display.
                Options: title, method, release_date, source, resolution
                Use 'none' to hide all header fields.
        input_form: If True, redirects to the input form with fields pre-filled
        base_url: Base URL for the R3DCID service
        
    Returns:
        Complete URL string with encoded query parameters.
        
    Examples:
        >>> build_r3dcid_url("2IZ8")
        'https://rna.bgsu.edu/fr3d/r3dcid?chains=2IZ8&format=html'
        
        >>> build_r3dcid_url(
        ...     "8GLP|1|L8+8GLP|1|L5",
        ...     coloring="grayscale",
        ...     dim=["bph", "br", "sr", "so"],
        ...     hide=["stacking", "near"],
        ...     text=["basepair", "bph", "br", "sr", "so"],
        ...     input_form=True
        ... )
        'https://rna.bgsu.edu/fr3d/r3dcid?chains=8GLP|1|L8+8GLP|1|L5&coloring=grayscale&dim=bph,br,sr,so&hide=stacking,near&text=basepair,bph,br,sr,so&input_form=True'
    """
    from urllib.parse import urlencode, quote
    
    # Helper to convert list or string to comma-separated string
    def to_csv(value: Optional[list[str] | str]) -> Optional[str]:
        if value is None:
            return None
        if isinstance(value, list):
            return ",".join(value) if value else None
        return value if value else None
    
    # Build params dict
    params: dict[str, str] = {"chains": chains}
    
    if assemblies:
        params["assemblies"] = assemblies
    
    if symmetries:
        params["symmetries"] = symmetries
    
    # Handle output format
    fmt = output_format.value if isinstance(output_format, OutputFormat) else output_format.lower()
    if fmt and fmt != "pdf":  # pdf is default, only include if different
        params["format"] = fmt
    
    # Handle coloring
    color = coloring.value if isinstance(coloring, ColorScheme) else coloring.lower()
    if color and color != "default":
        params["coloring"] = color
    
    if description:
        params["description"] = description
    
    # Handle show/dim/hide
    show_csv = to_csv(show)
    dim_csv = to_csv(dim)
    hide_csv = to_csv(hide)
    
    if show_csv:
        params["show"] = show_csv
    if dim_csv:
        params["dim"] = dim_csv
    if hide_csv:
        params["hide"] = hide_csv
    
    # Handle text options
    text_csv = to_csv(text)
    if text_csv and text_csv != "basepair":
        params["text"] = text_csv
    
    # Handle helix size
    if helix_size is not None:
        params["hs"] = str(helix_size)
    
    # Handle n3d
    if not n3d:
        params["n3d"] = "false"
    
    # Handle header fields
    header_csv = to_csv(header)
    if header_csv:
        default_headers = "title,method,release_date,source,resolution"
        if header_csv.lower() == "none" or header_csv == "":
            params["header"] = "none"
        elif sorted(header_csv.split(",")) != sorted(default_headers.split(",")):
            params["header"] = header_csv
    
    # Handle input_form
    if input_form:
        params["input_form"] = "True"
    
    # Build URL - use custom encoding to preserve | and + and , in chains
    query_parts = []
    for key, value in params.items():
        if key == "chains":
            # Don't encode | and + and , in chains parameter
            encoded_value = quote(value, safe="|+,")
        else:
            encoded_value = quote(str(value), safe=",")
        query_parts.append(f"{key}={encoded_value}")
    
    return f"{base_url}?{'&'.join(query_parts)}"


def build_r3dcid_url_from_model(
    input_data: R3DCIDInput,
    *,
    input_form: bool = False,
    base_url: str = DEFAULT_BASE_URL,
) -> str:
    """
    Build an R3DCID URL from a validated R3DCIDInput model.
    
    Args:
        input_data: Validated R3DCIDInput model instance.
        input_form: If True, redirects to the input form with fields pre-filled.
        base_url: Base URL for the R3DCID service.
        
    Returns:
        Complete URL string with encoded query parameters.
    """
    from urllib.parse import quote
    
    params = input_data.to_url_params()
    
    if input_form:
        params["input_form"] = "True"
    
    # Build URL
    query_parts = []
    for key, value in params.items():
        if key == "chains":
            encoded_value = quote(value, safe="|+,")
        else:
            encoded_value = quote(str(value), safe=",")
        query_parts.append(f"{key}={encoded_value}")
    
    return f"{base_url}?{'&'.join(query_parts)}"


def parse_r3dcid_url(url: str) -> R3DCIDQueryParams:
    """
    Parse an R3DCID URL and return a query params model.
    
    Args:
        url: Full R3DCID URL string.
        
    Returns:
        R3DCIDQueryParams model with parsed parameters.
        
    Example:
        >>> params = parse_r3dcid_url(
        ...     "https://rna.bgsu.edu/fr3d/r3dcid?chains=2IZ8&dim=nested-non-wc"
        ... )
        >>> params.chains
        '2IZ8'
        >>> params.dim
        'nested-non-wc'
    """
    from urllib.parse import urlparse, parse_qs, unquote
    
    parsed = urlparse(url)
    query_dict = parse_qs(parsed.query)
    
    # Convert from dict[str, list[str]] to dict[str, str]
    flat_params: dict[str, Optional[str]] = {}
    for key, values in query_dict.items():
        flat_params[key] = unquote(values[0]) if values else None
    
    # Ensure chains is present
    if "chains" not in flat_params or not flat_params["chains"]:
        raise ValueError("URL must contain a 'chains' parameter")
    
    return R3DCIDQueryParams(**flat_params)


# ============================================================================
# Main Entry Point (Examples)
# ============================================================================

if __name__ == "__main__":
    """
    Example usage demonstrating URL building and parsing.
    """
    
    print("=" * 70)
    print("R3DCID Abstract Module - Examples")
    print("=" * 70)
    
    # Example 1: Build URL using function with keyword arguments
    print("\n1. Building URL with build_r3dcid_url():")
    print("-" * 50)
    
    url1 = build_r3dcid_url(
        "7JQQ|1|K,L,M,N,O",
        coloring="wong",
        text=["basepair", "bph", "br", "so"],
        input_form=True
    )
    print(f"   URL: {url1}")
    
    # Example 2: Parse the URL back
    print("\n2. Parsing the URL with parse_r3dcid_url():")
    print("-" * 50)
    
    example_url = "https://rna.bgsu.edu/fr3d/r3dcid?chains=8GLP|1|L8+8GLP|1|L5&coloring=grayscale&dim=bph,br,sr,so&hide=stacking,near&text=basepair,bph,br,sr,so&input_form=True"
    params = parse_r3dcid_url(example_url)
    print(f"   chains: {params.chains}")
    print(f"   assemblies: {params.assemblies}")
    print(f"   symmetries: {params.symmetries}")
    print(f"   format: {params.format}")
    print(f"   coloring: {params.coloring}")
    print(f"   description: {params.description}")
    print(f"   show: {params.show}")
    print(f"   dim: {params.dim}")
    print(f"   hide: {params.hide}")
    print(f"   text: {params.text}")
    print(f"   hs: {params.hs}")
    print(f"   n3d: {params.n3d}")
    print(f"   header: {params.header}")
    print(f"   input_form: {params.input_form}")
    
    # Example 3: Convert to full R3DCIDInput model
    print("\n3. Converting to R3DCIDInput model:")
    print("-" * 50)
    
    input_model = params.to_r3dcid_input()
    print(f"   chains: {input_model.chains}")
    print(f"   assemblies: {input_model.assemblies}")
    print(f"   symmetries: {input_model.symmetries}")
    print(f"   output_format: {input_model.output_format}")
    print(f"   color_scheme: {input_model.color_scheme}")
    print(f"   description: {input_model.description}")
    print(f"   text_options: {input_model.text_options}")
    print(f"   helix_size: {input_model.helix_size}")
    print(f"   show_no_3d_coords: {input_model.show_no_3d_coords}")
    print(f"   header_fields: {input_model.header_fields}")
    print(f"   interaction_arcs:")
    print(f"      nested_wc: {input_model.interaction_arcs.nested_wc}")
    print(f"      lr_wc: {input_model.interaction_arcs.lr_wc}")
    print(f"      nested_non_wc: {input_model.interaction_arcs.nested_non_wc}")
    print(f"      lr_non_wc: {input_model.interaction_arcs.lr_non_wc}")
    print(f"      stacking: {input_model.interaction_arcs.stacking}")
    print(f"      bph: {input_model.interaction_arcs.bph}")
    print(f"      br: {input_model.interaction_arcs.br}")
    print(f"      sr: {input_model.interaction_arcs.sr}")
    print(f"      so: {input_model.interaction_arcs.so}")
    print(f"      near: {input_model.interaction_arcs.near}")
    
    # Example 4: Build URL from validated model
    print("\n4. Building URL from R3DCIDInput model:")
    print("-" * 50)
    
    url2 = build_r3dcid_url_from_model(input_model, input_form=True)
    print(f"   URL: {url2}")
    
    # Example 5: Create model directly and build URL
    print("\n5. Creating R3DCIDInput directly:")
    print("-" * 50)
    
    custom_input = R3DCIDInput(
        chains="8GLP|1|L8+8GLP|1|L5",
        color_scheme=ColorScheme.GRAYSCALE,
        text_options=[TextOption.BASEPAIR, TextOption.BPH, TextOption.BR, TextOption.SR, TextOption.SO],
        show_no_3d_coords=False,
    )
    custom_input.interaction_arcs.stacking = InteractionVisibility.HIDE
    custom_input.interaction_arcs.near = InteractionVisibility.HIDE
    custom_input.interaction_arcs.bph = InteractionVisibility.DIM
    custom_input.interaction_arcs.br = InteractionVisibility.DIM
    custom_input.interaction_arcs.sr = InteractionVisibility.DIM
    custom_input.interaction_arcs.so = InteractionVisibility.DIM
    
    url3 = build_r3dcid_url_from_model(custom_input, input_form=True)
    print(f"   URL: {url3}")
    
    # Example 6: Simple PDB ID only
    print("\n6. Simple PDB ID only:")
    print("-" * 50)
    
    url4 = build_r3dcid_url("7EZ2", helix_size=9)
    print(f"   URL: {url4}")
    
    # Example 7: With assemblies
    print("\n7. With assemblies:")
    print("-" * 50)
    
    url5 = build_r3dcid_url("4V9O", assemblies="4,2")
    print(f"   URL: {url5}")
    
    print("\n" + "=" * 70)
    print("Examples completed successfully!")
    print("=" * 70)
