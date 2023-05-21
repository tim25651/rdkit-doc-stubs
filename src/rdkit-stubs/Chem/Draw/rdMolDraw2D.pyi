"""
rdkit.Chem.Draw.rdMolDraw2D module¶
Module containing a C++ implementation of 2D molecule drawing
"""
from typing import Any, ClassVar, overload

import Boost.Python
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdChemReactions import ChemicalReaction
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters
from rdkit.Geometry.rdGeometry import Point2D

class ContourParams(Boost.Python.instance):
    """
    Parameters for drawing contours

    C++ signature :void __init__(_object*)

    property contourWidth¶
    line width of the contours

    property dashNegative¶
    use a dashed line for negative contours

    property extraGridPadding¶
    extra space (in molecule coords) around the grid

    property fillGrid¶
    colors the grid in addition to drawing contours

    property gridResolution¶
    set the resolution of the grid"""

    __instance_size__: ClassVar[int] = ...
    contourWidth: Any
    dashNegative: Any
    extraGridPadding: Any
    fillGrid: Any
    gridResolution: Any
    setScale: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def setColourMap(self: ContourParams, colours: AtomPairsParameters) -> None:
        """
        C++ signature :void setColourMap(RDKit::MolDraw2DUtils::ContourParams {lvalue},boost::python::api::object)
        """
        ...
    def setContourColour(self: ContourParams, colour: tuple) -> None:
        """
        C++ signature :void setContourColour(RDKit::MolDraw2DUtils::ContourParams {lvalue},boost::python::tuple)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class IntStringMap(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    ...
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __contains__(cls, other) -> Any: ...
    @classmethod
    def __delitem__(cls, other) -> Any: ...
    @classmethod
    def __getitem__(cls, index) -> Any: ...
    @classmethod
    def __iter__(cls) -> Any: ...
    @classmethod
    def __len__(cls) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, index, object) -> Any: ...

class MolDraw2D(Boost.Python.instance):
    """
    Drawer abstract base class
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def ClearDrawing(self: MolDraw2D) -> None:
        """
        clears the drawing by filling it with the background color

        C++ signature :void ClearDrawing(RDKit::MolDraw2D {lvalue})"""
        ...
    def DrawArc(
        self: MolDraw2D,
        center: Point2D,
        radius: float,
        angle1: float,
        angle2: float,
        rawCoords: bool = False,
    ) -> None:
        """
        draws an arc with the current drawing style. The coordinates are in the molecule frame, the angles are in degrees, angle2 should be > angle1.

        C++ signature :void DrawArc(RDKit::MolDraw2D {lvalue},RDGeom::Point2D,double,double,double [,bool=False])
        """
        ...
    def DrawArrow(
        self: MolDraw2D,
        cds1: Point2D,
        cds2: Point2D,
        asPolygon: bool = False,
        frac: float = 0.05,
        angle: float = 0.5235987755982988,
        color: AtomPairsParameters = None,
        rawCoords: bool = False,
    ) -> None:
        """
        draws an arrow with the current drawing style. The coordinates are in the molecule frame. If asPolygon is true the head of the arrow will be drawn as a triangle, otherwise two lines are used.

        C++ signature :void DrawArrow(RDKit::MolDraw2D {lvalue},RDGeom::Point2D,RDGeom::Point2D [,bool=False [,double=0.05 [,double=0.5235987755982988 [,boost::python::api::object=None [,bool=False]]]]])
        """
        ...
    def DrawAttachmentLine(
        self: MolDraw2D,
        cds1: Point2D,
        cds2: Point2D,
        color: tuple,
        len: float = 1.0,
        nSegments: int = 16,
        rawCoords: bool = False,
    ) -> None:
        """
        draw a line indicating the presence of an attachment point (normally a squiggle line perpendicular to a bond)

        C++ signature :void DrawAttachmentLine(RDKit::MolDraw2D {lvalue},RDGeom::Point2D,RDGeom::Point2D,boost::python::tuple {lvalue} [,double=1.0 [,unsigned int=16 [,bool=False]]])
        """
        ...
    def DrawEllipse(
        self: MolDraw2D, cds1: Point2D, cds2: Point2D, rawCoords: bool = False
    ) -> None:
        """
        draws a triangle with the current drawing style in the rectangle defined by the two points. The coordinates are in the molecule frame

        C++ signature :void DrawEllipse(RDKit::MolDraw2D {lvalue},RDGeom::Point2D,RDGeom::Point2D [,bool=False])
        """
        ...
    def DrawLine(
        self: MolDraw2D, cds1: Point2D, cds2: Point2D, rawCoords: bool = False
    ) -> None:
        """
        draws a line with the current drawing style. The coordinates are in the molecule frame

        C++ signature :void DrawLine(RDKit::MolDraw2D {lvalue},RDGeom::Point2D,RDGeom::Point2D [,bool=False])
        """
        ...
    @overload
    def DrawMolecule(
        self: MolDraw2D,
        mol: Mol,
        highlightAtoms: AtomPairsParameters = None,
        highlightAtomColors: AtomPairsParameters = None,
        highlightAtomRadii: AtomPairsParameters = None,
        confId: int = -1,
        legend: str = "",
    ) -> None:
        """
        renders a molecule

        C++ signature :void DrawMolecule(RDKit::MolDraw2D {lvalue},RDKit::ROMol,boost::python::api::object,boost::python::api::object [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,int=-1 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’]]]]])
        """
        ...
    @overload
    def DrawMolecule(
        self: MolDraw2D,
        mol: Mol,
        highlightAtoms: AtomPairsParameters,
        highlightBonds: AtomPairsParameters,
        highlightAtomColors: AtomPairsParameters = None,
        highlightBondColors: AtomPairsParameters = None,
        highlightAtomRadii: AtomPairsParameters = None,
        confId: int = -1,
        legend: str = "",
    ) -> None: ...
    def DrawMoleculeWithHighlights(
        self: MolDraw2D,
        mol: Mol,
        legend: str,
        highlight_atom_map: AtomPairsParameters,
        highlight_bond_map: AtomPairsParameters,
        highlight_radii: AtomPairsParameters,
        highlight_linewidth_multipliers: AtomPairsParameters,
        confId: int = -1,
    ) -> None:
        """
        renders a molecule with multiple highlight colours

        C++ signature :void DrawMoleculeWithHighlights(RDKit::MolDraw2D {lvalue},RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object [,int=-1])
        """
        ...
    def DrawMolecules(
        self: MolDraw2D,
        mols: AtomPairsParameters,
        highlightAtoms: AtomPairsParameters = None,
        highlightBonds: AtomPairsParameters = None,
        highlightAtomColors: AtomPairsParameters = None,
        highlightBondColors: AtomPairsParameters = None,
        highlightAtomRadii: AtomPairsParameters = None,
        confIds: AtomPairsParameters = None,
        legends: AtomPairsParameters = None,
    ) -> None:
        """
        renders multiple molecules

        C++ signature :void DrawMolecules(RDKit::MolDraw2D {lvalue},boost::python::api::object [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None]]]]]]])
        """
        ...
    def DrawPolygon(
        self: MolDraw2D, cds: AtomPairsParameters, rawCoords: bool = False
    ) -> None:
        """
        draws a polygon with the current drawing style. The coordinates are in the molecule frame

        C++ signature :void DrawPolygon(RDKit::MolDraw2D {lvalue},boost::python::api::object [,bool=False])
        """
        ...
    def DrawReaction(
        self: MolDraw2D,
        rxn: ChemicalReaction,
        highlightByReactant: bool = False,
        highlightColorsReactants: AtomPairsParameters = None,
        confIds: AtomPairsParameters = None,
    ) -> None:
        """
        renders a reaction

        C++ signature :void DrawReaction(RDKit::MolDraw2D {lvalue},RDKit::ChemicalReaction [,bool=False [,boost::python::api::object=None [,boost::python::api::object=None]]])
        """
        ...
    def DrawRect(
        self: MolDraw2D, cds1: Point2D, cds2: Point2D, rawCoords: bool = False
    ) -> None:
        """
        draws a rectangle with the current drawing style in the rectangle defined by the two points. The coordinates are in the molecule frame

        C++ signature :void DrawRect(RDKit::MolDraw2D {lvalue},RDGeom::Point2D,RDGeom::Point2D [,bool=False])
        """
        ...
    @overload
    def DrawString(
        self: MolDraw2D, string: str, pos: Point2D, rawCoords: bool = False
    ) -> None:
        """
        add aligned text to the canvas. The align argument can be 0 (=MIDDLE), 1 (=START), or 2 (=END)

        C++ signature :void DrawString(RDKit::MolDraw2D {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDGeom::Point2D,int [,bool=False])
        """
        ...
    @overload
    def DrawString(
        self: MolDraw2D, string: str, pos: Point2D, align: int, rawCoords: bool = False
    ) -> None: ...
    def DrawTriangle(
        self: MolDraw2D,
        cds1: Point2D,
        cds2: Point2D,
        cds3: Point2D,
        rawCoords: bool = False,
    ) -> None:
        """
        draws a triangle with the current drawing style. The coordinates are in the molecule frame

        C++ signature :void DrawTriangle(RDKit::MolDraw2D {lvalue},RDGeom::Point2D,RDGeom::Point2D,RDGeom::Point2D [,bool=False])
        """
        ...
    def DrawWavyLine(
        self: MolDraw2D,
        cds1: Point2D,
        cds2: Point2D,
        color1: tuple,
        color2: tuple,
        nSegments: int = 16,
        vertOffset: float = 0.05,
        rawCoords: bool = False,
    ) -> None:
        """
        draw a line indicating the presence of an attachment point (normally a squiggle line perpendicular to a bond)

        C++ signature :void DrawWavyLine(RDKit::MolDraw2D {lvalue},RDGeom::Point2D,RDGeom::Point2D,boost::python::tuple {lvalue},boost::python::tuple {lvalue} [,unsigned int=16 [,double=0.05 [,bool=False]]])
        """
        ...
    def FillPolys(self, arg1: MolDraw2D) -> bool:
        """
        returns whether or not polygons are being filled

        C++ signature :bool FillPolys(RDKit::MolDraw2D {lvalue})"""
        ...
    def FlexiMode(self, arg1: MolDraw2D) -> bool:
        """
        returns whether or not FlexiMode is being used

        C++ signature :bool FlexiMode(RDKit::MolDraw2D {lvalue})"""
        ...
    def FontSize(self, arg1: MolDraw2D) -> float:
        """
        get the default font size. The units are, roughly, pixels.

        C++ signature :double FontSize(RDKit::MolDraw2D {lvalue})"""
        ...
    @overload
    def GetDrawCoords(self: MolDraw2D, point: Point2D) -> Point2D:
        """
        get the coordinates in drawing space for a particular atom

        C++ signature :RDGeom::Point2D GetDrawCoords(RDKit::MolDraw2D {lvalue},int)"""
        ...
    @overload
    def GetDrawCoords(self: MolDraw2D, atomIndex: int) -> Point2D: ...
    def GetMolSize(
        self: MolDraw2D,
        mol: Mol,
        highlightAtoms: AtomPairsParameters = None,
        highlightBonds: AtomPairsParameters = None,
        highlightAtomColors: AtomPairsParameters = None,
        highlightBondColors: AtomPairsParameters = None,
        highlightAtomRadii: AtomPairsParameters = None,
        confId: int = -1,
        legend: str = "",
    ) -> tuple:
        """
        returns the width and height required to draw a molecule at the current size

        C++ signature :boost::python::tuple GetMolSize(RDKit::MolDraw2D {lvalue},RDKit::ROMol [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,int=-1 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’]]]]]]])
        """
        ...
    def Height(self, arg1: MolDraw2D) -> int:
        """
        get the height of the drawing canvas

        C++ signature :int Height(RDKit::MolDraw2D {lvalue})"""
        ...
    def LineWidth(self, arg1: MolDraw2D) -> float:
        """
        returns the line width being used

        C++ signature :double LineWidth(RDKit::MolDraw2D {lvalue})"""
        ...
    def Offset(self, arg1: MolDraw2D) -> Point2D:
        """
        returns the offset (in drawing coordinates) for the drawing

        C++ signature :RDGeom::Point2D Offset(RDKit::MolDraw2D {lvalue})"""
        ...
    def SetColour(self, arg1: MolDraw2D, arg2: tuple) -> None:
        """
        set the color being used fr drawing and filling

        C++ signature :void SetColour(RDKit::MolDraw2D {lvalue},boost::python::tuple)"""
        ...
    def SetDrawOptions(self, arg1: MolDraw2D, arg2: MolDrawOptions) -> None:
        """
        Copies the drawing options passed in over our drawing options

        C++ signature :void SetDrawOptions(RDKit::MolDraw2D {lvalue},RDKit::MolDrawOptions)
        """
        ...
    def SetFillPolys(self, arg1: MolDraw2D, arg2: bool) -> None:
        """
        sets whether or not polygons are filled

        C++ signature :void SetFillPolys(RDKit::MolDraw2D {lvalue},bool)"""
        ...
    def SetFlexiMode(self, arg1: MolDraw2D, arg2: bool) -> None:
        """
        when FlexiMode is set, molecules will always been drawn with the default values for bond length, font size, etc.

        C++ signature :void SetFlexiMode(RDKit::MolDraw2D {lvalue},bool)"""
        ...
    def SetFontSize(self, arg1: MolDraw2D, arg2: float) -> None:
        """
        change the default font size. The units are, roughly, pixels.

        C++ signature :void SetFontSize(RDKit::MolDraw2D {lvalue},double)"""
        ...
    def SetLineWidth(self, arg1: MolDraw2D, arg2: float) -> None:
        """
        set the line width being used

        C++ signature :void SetLineWidth(RDKit::MolDraw2D {lvalue},double)"""
        ...
    def SetOffset(self, arg1: MolDraw2D, arg2: int, arg3: int) -> None:
        """
        set the offset (in drawing coordinates) for the drawing

        C++ signature :void SetOffset(RDKit::MolDraw2D {lvalue},int,int)"""
        ...
    def SetScale(
        self: MolDraw2D,
        width: int,
        height: int,
        minv: Point2D,
        maxv: Point2D,
        mol: AtomPairsParameters = None,
    ) -> None:
        """
        uses the values provided to set the drawing scaling

        C++ signature :void SetScale(RDKit::MolDraw2D {lvalue},int,int,RDGeom::Point2D,RDGeom::Point2D [,boost::python::api::object=None])
        """
        ...
    def Width(self, arg1: MolDraw2D) -> int:
        """
        get the width of the drawing canvas

        C++ signature :int Width(RDKit::MolDraw2D {lvalue})"""
        ...
    def drawOptions(self, arg1: MolDraw2D) -> MolDrawOptions:
        """
        Returns a modifiable version of the current drawing options

        C++ signature :RDKit::MolDrawOptions {lvalue} drawOptions(RDKit::MolDraw2D {lvalue})
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MolDraw2DCairo(MolDraw2D):
    """
    Cairo molecule drawer

    C++ signature :void __init__(_object*,int,int [,int=-1 [,int=-1 [,bool=False]]])"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def FinishDrawing(self, arg1: MolDraw2DCairo) -> None:
        """
        add the last bits to finish the drawing

        C++ signature :void FinishDrawing(RDKit::MolDraw2DCairo {lvalue})"""
        ...
    def GetDrawingText(self, arg1: MolDraw2DCairo) -> object:
        """
        return the PNG data as a string

        C++ signature :boost::python::api::object GetDrawingText(RDKit::MolDraw2DCairo)
        """
        ...
    def WriteDrawingText(self, arg1: MolDraw2DCairo, arg2: str) -> None:
        """
        write the PNG data to the named file

        C++ signature :void WriteDrawingText(RDKit::MolDraw2DCairo {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MolDraw2DSVG(MolDraw2D):
    """
    SVG molecule drawer

    C++ signature :void __init__(_object*,int,int [,int=-1 [,int=-1 [,bool=False]]])"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AddMoleculeMetadata(
        self, arg1: MolDraw2DSVG, mol: Mol, confId: int = -1
    ) -> None:
        """
        add RDKit-specific information to the bottom of the drawing

        C++ signature :void AddMoleculeMetadata(RDKit::MolDraw2DSVG {lvalue},RDKit::ROMol [,int=-1])
        """
        ...
    def FinishDrawing(self, arg1: MolDraw2DSVG) -> None:
        """
        add the last bits of SVG to finish the drawing

        C++ signature :void FinishDrawing(RDKit::MolDraw2DSVG {lvalue})"""
        ...
    def GetDrawingText(self, arg1: MolDraw2DSVG) -> str:
        """
        return the SVG

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetDrawingText(RDKit::MolDraw2DSVG {lvalue})
        """
        ...
    def TagAtoms(
        self,
        arg1: MolDraw2DSVG,
        mol: Mol,
        radius: float = 0.2,
        events: AtomPairsParameters = None,
    ) -> None:
        """
        allow atom selection in the SVG

        C++ signature :void TagAtoms(RDKit::MolDraw2DSVG {lvalue},RDKit::ROMol [,double=0.2 [,boost::python::api::object=None]])
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MolDrawOptions(Boost.Python.instance):
    """
    Drawing options

    C++ signature :void __init__(_object*)

    property addAtomIndices¶
    adds atom indices to drawings. Default False.

    property addBondIndices¶
    adds bond indices to drawings. Default False.

    property addStereoAnnotation¶
    adds R/S and E/Z to drawings. Default False.

    property additionalAtomLabelPadding¶
    additional padding to leave around atom labels. Expressed as a fraction of the font size.

    property annotationFontScale¶
    Scale of font for atom and bond annotation relative to atomlabel font.  Default=0.75.

    property atomHighlightsAreCircles¶
    forces atom highlights always to be circles.Default (false) is to put ellipses roundlonger labels.

    property atomLabelDeuteriumTritium¶
    labels deuterium as D and tritium as T

    property atomLabels¶
    maps indices to atom labels

    property atomRegions¶
    regions to outline

    property baseFontSize¶
    relative size of font.  Defaults to 0.6.  -1 means use default.

    property bondLineWidth¶
    if positive, this overrides the default line width for bonds

    property centreMoleculesBeforeDrawing¶
    Moves the centre of the drawn molecule to (0,0).Default False.

    property circleAtoms¶

    property clearBackground¶
    clear the background before drawing a molecule

    property comicMode¶
    simulate hand-drawn lines for bonds. When combined with a font like Comic-Sans or Comic-Neue, this gives xkcd-like drawings. Default is false.

    property continuousHighlight¶

    property drawMolsSameScale¶
    when drawing multiple molecules with DrawMolecules, forces them to use the same scale.  Default is true.

    property dummiesAreAttachments¶

    property dummyIsotopeLabels¶
    adds isotope labels on dummy atoms. Default True.

    property explicitMethyl¶
    Draw terminal methyls explictly.  Default is false.

    property fillHighlights¶

    property fixedBondLength¶
    If > 0.0, fixes bond length to this number of pixelsunless that would make it too big.  Default -1.0 meansno fix.  If both set, fixedScale takes precedence.

    property fixedFontSize¶
    font size in pixels. default=-1 means not fixed.  If set, always used irrespective of scale, minFontSize and maxFontSize.

    property fixedScale¶
    If > 0.0, fixes scale to that fraction of width ofdraw window.  Default -1.0 means adjust scale to fit.

    property flagCloseContactsDist¶

    property fontFile¶
    Font file for use with FreeType text drawer.  Can also be BuiltinTelexRegular (the default) or BuiltinRobotoRegular.
    """

    __instance_size__: ClassVar[int] = ...
    addAtomIndices: Any
    addBondIndices: Any
    addStereoAnnotation: Any
    additionalAtomLabelPadding: Any
    annotationFontScale: Any
    atomHighlightsAreCircles: Any
    atomLabelDeuteriumTritium: Any
    atomLabels: Any
    atomRegions: Any
    baseFontSize: Any
    bondLineWidth: Any
    centreMoleculesBeforeDrawing: Any
    circleAtoms: Any
    clearBackground: Any
    comicMode: Any
    continuousHighlight: Any
    drawMolsSameScale: Any
    dummiesAreAttachments: Any
    dummyIsotopeLabels: Any
    explicitMethyl: Any
    fillHighlights: Any
    fixedBondLength: Any
    fixedFontSize: Any
    fixedScale: Any
    flagCloseContactsDist: Any
    fontFile: Any
    highlightBondWidthMultiplier: Any
    highlightRadius: Any
    includeAtomTags: Any
    includeChiralFlagLabel: Any
    includeMetadata: Any
    includeRadicals: Any
    isotopeLabels: Any
    legendFontSize: Any
    legendFraction: Any
    maxFontSize: Any
    minFontSize: Any
    multipleBondOffset: Any
    noAtomLabels: Any
    padding: Any
    prepareMolsBeforeDrawing: Any
    rotate: Any
    scaleBondWidth: Any
    scaleHighlightBondWidth: Any
    scalingFactor: Any
    simplifiedStereoGroupLabel: Any
    singleColourWedgeBonds: Any
    splitBonds: Any
    unspecifiedStereoIsUnknown: Any
    useComplexQueryAtomSymbols: Any
    useMolBlockWedging: Any
    variableAtomRadius: Any
    variableBondWidthMultiplier: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def getAnnotationColour(self, arg1: MolDrawOptions) -> object:
        """
        method returning the annotation colour

        C++ signature :boost::python::api::object getAnnotationColour(RDKit::MolDrawOptions)
        """
        ...
    def getBackgroundColour(self, arg1: MolDrawOptions) -> object:
        """
        method returning the background colour

        C++ signature :boost::python::api::object getBackgroundColour(RDKit::MolDrawOptions)
        """
        ...
    def getHighlightColour(self, arg1: MolDrawOptions) -> object:
        """
        method returning the highlight colour

        C++ signature :boost::python::api::object getHighlightColour(RDKit::MolDrawOptions)
        """
        ...
    def getLegendColour(self, arg1: MolDrawOptions) -> object:
        """
        method returning the legend colour

        C++ signature :boost::python::api::object getLegendColour(RDKit::MolDrawOptions)
        """
        ...
    def getQueryColour(self, arg1: MolDrawOptions) -> object:
        """
        method returning the query colour

        C++ signature :boost::python::api::object getQueryColour(RDKit::MolDrawOptions)
        """
        ...
    def getSymbolColour(self, arg1: MolDrawOptions) -> object:
        """
        method returning the symbol colour

        C++ signature :boost::python::api::object getSymbolColour(RDKit::MolDrawOptions)
        """
        ...
    def getVariableAttachmentColour(self, arg1: MolDrawOptions) -> object:
        """
        method for getting the colour of variable attachment points

        C++ signature :boost::python::api::object getVariableAttachmentColour(RDKit::MolDrawOptions)
        """
        ...
    def setAnnotationColour(self, arg1: MolDrawOptions, arg2: tuple) -> None:
        """
        method for setting the annotation colour

        C++ signature :void setAnnotationColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)
        """
        ...
    def setAtomPalette(self, arg1: MolDrawOptions, arg2: AtomPairsParameters) -> None:
        """
        sets the palette for atoms and bonds from a dictionary mapping ints to 3-tuples

        C++ signature :void setAtomPalette(RDKit::MolDrawOptions {lvalue},boost::python::api::object)
        """
        ...
    def setBackgroundColour(self, arg1: MolDrawOptions, arg2: tuple) -> None:
        """
        method for setting the background colour

        C++ signature :void setBackgroundColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)
        """
        ...
    def setHighlightColour(self, arg1: MolDrawOptions, arg2: tuple) -> None:
        """
        method for setting the highlight colour

        C++ signature :void setHighlightColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)
        """
        ...
    def setLegendColour(self, arg1: MolDrawOptions, arg2: tuple) -> None:
        """
        method for setting the legend colour

        C++ signature :void setLegendColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)
        """
        ...
    def setQueryColour(self, arg1: MolDrawOptions, arg2: tuple) -> None:
        """
        method for setting the query colour

        C++ signature :void setQueryColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)
        """
        ...
    def setSymbolColour(self, arg1: MolDrawOptions, arg2: tuple) -> None:
        """
        method for setting the symbol colour

        C++ signature :void setSymbolColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)
        """
        ...
    def setVariableAttachmentColour(self, arg1: MolDrawOptions, arg2: tuple) -> None:
        """
        method for setting the colour of variable attachment points

        C++ signature :void setVariableAttachmentColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)
        """
        ...
    def updateAtomPalette(
        self, arg1: MolDrawOptions, arg2: AtomPairsParameters
    ) -> None:
        """
        updates the palette for atoms and bonds from a dictionary mapping ints to 3-tuples

        C++ signature :void updateAtomPalette(RDKit::MolDrawOptions {lvalue},boost::python::api::object)
        """
        ...
    def useAvalonAtomPalette(self, arg1: MolDrawOptions) -> None:
        """
        use the Avalon renderer palette for atoms and bonds

        C++ signature :void useAvalonAtomPalette(RDKit::MolDrawOptions {lvalue})"""
        ...
    def useBWAtomPalette(self, arg1: MolDrawOptions) -> None:
        """
        use a black and white palette for atoms and bonds

        C++ signature :void useBWAtomPalette(RDKit::MolDrawOptions {lvalue})"""
        ...
    def useCDKAtomPalette(self, arg1: MolDrawOptions) -> None:
        """
        use the CDK palette for atoms and bonds

        C++ signature :void useCDKAtomPalette(RDKit::MolDrawOptions {lvalue})"""
        ...
    def useDefaultAtomPalette(self, arg1: MolDrawOptions) -> None:
        """
        use the default colour palette for atoms and bonds

        C++ signature :void useDefaultAtomPalette(RDKit::MolDrawOptions {lvalue})"""
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class map_indexing_suite_IntStringMap_entry(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def data(self, arg1: map_indexing_suite_IntStringMap_entry) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > data(std::pair<int const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > {lvalue})
        """
        ...
    def key(self, arg1: map_indexing_suite_IntStringMap_entry) -> int:
        """
        C++ signature :int key(std::pair<int const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > {lvalue})
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

def ContourAndDrawGaussians(
    self,
    drawer: MolDraw2D,
    locs: AtomPairsParameters,
    heights: AtomPairsParameters,
    widths: AtomPairsParameters,
    nContours: int = 10,
    levels: AtomPairsParameters = None,
    params: ContourParams = ...,
    mol: AtomPairsParameters = None,
) -> None:
    """
    Generates and draws contours for a set of gaussians

    drawer: the MolDraw2D object to use
    locs: locations of the gaussians
    heights: the heights (or weights) of the gaussians
    widths: the standard deviations of the gaussians
    nContours: the number of contours to draw
    levels: the contours to use
    ps: additional parameters controlling the contouring.
    mol: molecule used to help set scale.

    The values are calculated on a grid with spacing params.gridResolution.
    If params.setScale  is set, the grid size will be calculated based on the
    locations of the gaussians and params.extraGridPadding. Otherwise the current
    size of the viewport will be used.
    If the levels argument is empty, the contour levels will be determined
    automatically from the max and min values on the grid and levels will
    be updated to include the contour levels.
    If params.fillGrid is set, the data on the grid will also be drawn using
    the color scheme in params.colourMap
    If mol is not 0, uses the molecule to help set the scale, assuming that
    it will be drawn over the plot, so needs to fit on it.

    */

    C++ signature :void ContourAndDrawGaussians(RDKit::MolDraw2D {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object [,unsigned int=10 [,boost::python::api::object=None [,RDKit::MolDraw2DUtils::ContourParams=<rdkit.Chem.Draw.rdMolDraw2D.ContourParams object at 0x7fd4801f2ac0> [,boost::python::api::object=None]]]])
    """
    ...

def ContourAndDrawGrid(
    self,
    drawer: MolDraw2D,
    data: AtomPairsParameters,
    xcoords: AtomPairsParameters,
    ycoords: AtomPairsParameters,
    nContours: int = 10,
    levels: AtomPairsParameters = None,
    params: ContourParams = ...,
    mol: AtomPairsParameters = None,
) -> None:
    """
    Generates and draws contours for data on a grid

    drawer: the MolDraw2D object to use
    data: numpy array with the data to be contoured
    xcoords: the x coordinates of the grid
    ycoords: the y coordinates of the grid
    nContours: the number of contours to draw
    levels: the contours to use
    ps: additional parameters controlling the contouring
    mol: molecule used to help set scale.

    The values are calculated on a grid with spacing params.gridResolution.
    If params.setScale  is set, the grid size will be calculated based on the
    locations of the gaussians and params.extraGridPadding. Otherwise the current
    size of the viewport will be used.
    If the levels argument is empty, the contour levels will be determined
    automatically from the max and min values on the grid and levels will
    be updated to include the contour levels.
    If params.fillGrid is set, the data on the grid will also be drawn using
    the color scheme in params.colourMap
    If mol is not 0, uses the molecule to help set the scale, assuming that
    it will be drawn over the plot, so needs to fit on it.

    */

    C++ signature :void ContourAndDrawGrid(RDKit::MolDraw2D {lvalue},boost::python::api::object {lvalue},boost::python::api::object {lvalue},boost::python::api::object {lvalue} [,unsigned int=10 [,boost::python::api::object {lvalue}=None [,RDKit::MolDraw2DUtils::ContourParams=<rdkit.Chem.Draw.rdMolDraw2D.ContourParams object at 0x7fd4801f2b80> [,boost::python::api::object=None]]]])
    """
    ...

def DrawMoleculeACS1996(
    self,
    drawer: MolDraw2D,
    mol: Mol,
    legend: str = "",
    highlightAtoms: AtomPairsParameters = None,
    highlightBonds: AtomPairsParameters = None,
    highlightAtomColors: AtomPairsParameters = None,
    highlightBondColors: AtomPairsParameters = None,
    highlightAtomRadii: AtomPairsParameters = None,
    confId: int = -1,
) -> None:
    """
    Draws molecule in ACS 1996 mode.

    C++ signature :void DrawMoleculeACS1996(RDKit::MolDraw2D {lvalue},RDKit::ROMol [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’ [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,int=-1]]]]]]])
    """
    ...

def MeanBondLength(self, mol: Mol, confId: int = -1) -> float:
    """
    Calculate the mean bond length for the molecule.

    C++ signature :double MeanBondLength(RDKit::ROMol [,int=-1])"""
    ...

def MolToACS1996SVG(
    self,
    mol: Mol,
    legend: str = "",
    highlightAtoms: AtomPairsParameters = None,
    highlightBonds: AtomPairsParameters = None,
    highlightAtomColors: AtomPairsParameters = None,
    highlightBondColors: AtomPairsParameters = None,
    highlightAtomRadii: AtomPairsParameters = None,
    confId: int = -1,
) -> str:
    """
    Returns ACS 1996 mode svg for a molecule

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToACS1996SVG(RDKit::ROMol [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’ [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,int=-1]]]]]]])
    """
    ...

def MolToSVG(
    self,
    mol: Mol,
    width: int = 300,
    height: int = 300,
    highlightAtoms: AtomPairsParameters = None,
    kekulize: bool = True,
    lineWidthMult: int = 1,
    fontSize: bool = 12,
    includeAtomCircles: int = True,
) -> str:
    """
    Returns svg for a molecule

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToSVG(RDKit::ROMol [,unsigned int=300 [,unsigned int=300 [,boost::python::api::object=None [,bool=True [,unsigned int=1 [,bool=12 [,int=True]]]]]]])
    """
    ...

def PrepareAndDrawMolecule(
    self,
    drawer: MolDraw2D,
    mol: Mol,
    legend: str = "",
    highlightAtoms: AtomPairsParameters = None,
    highlightBonds: AtomPairsParameters = None,
    highlightAtomColors: AtomPairsParameters = None,
    highlightBondColors: AtomPairsParameters = None,
    highlightAtomRadii: AtomPairsParameters = None,
    confId: int = -1,
    kekulize: bool = True,
) -> None:
    """
    Preps a molecule for drawing and actually draws it

    C++ signature :void PrepareAndDrawMolecule(RDKit::MolDraw2D {lvalue},RDKit::ROMol [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’ [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,int=-1 [,bool=True]]]]]]]])
    """
    ...

def PrepareMolForDrawing(
    self,
    mol: Mol,
    kekulize: bool = True,
    addChiralHs: bool = True,
    wedgeBonds: bool = True,
    forceCoords: bool = False,
    wavyBonds: bool = False,
) -> Mol:
    """
    Does some cleanup operations on the molecule to prepare it to draw nicely.
    The operations include: kekulization, addition of chiral Hs (so that we can draw
    wedges to them), wedging of bonds at chiral centers, and generation of a 2D
    conformation if the molecule does not already have a conformation
    Returns a modified copy of the molecule.

    C++ signature :RDKit::ROMol* PrepareMolForDrawing(RDKit::ROMol const* [,bool=True [,bool=True [,bool=True [,bool=False [,bool=False]]]]])
    """
    ...

def SetACS1996Mode(self, drawOptions: MolDrawOptions, meanBondLength: float) -> None:
    """
    Set the draw options to produce something as close as possible to
    the ACS 1996 guidelines as described at
    https://en.wikipedia.org/wiki/Wikipedia:Manual_of_Style/Chemistry/Structure_drawing

    MolDrawOptions opt - the options what will be changed
    float meanBondLength - mean bond length of the molecule

    Works best if the MolDraw2D object is created with width and height -1 (a
    flexiCanvas).
    The mean bond length may be calculated with MeanBondLength.
    It is used to calculate the offset for the lines in multiple bonds.

    Options changed are:bondLineWidth = 0.6
    scaleBondWidth = false
    scalingFactor = 14.4 / meanBondLen
    multipleBondOffset = 0.18
    highlightBondWidthMultiplier = 32
    setMonochromeMode - black and white
    fixedFontSize = 10
    additionalAtomLabelPadding = 0.066
    fontFile - if it isn’t set already, then if RDBASE is set and the file

    exists, uses $RDBASE/Data/Fonts/FreeSans.ttf.  Otherwise uses
    BuiltinRobotoRegular.

    */

    C++ signature :void SetACS1996Mode(RDKit::MolDrawOptions {lvalue},double)"""
    ...

@overload
def SetDarkMode(self, arg1: MolDrawOptions) -> None:
    """
    set dark mode for a MolDraw2D object

    C++ signature :void SetDarkMode(RDKit::MolDraw2D {lvalue})"""
    ...

@overload
def SetDarkMode(self, arg1: MolDraw2D) -> None: ...
@overload
def SetMonochromeMode(
    self, options: MolDrawOptions, fgColour: tuple, bgColour: tuple
) -> None:
    """
    set monochrome mode for a MolDraw2D object

    C++ signature :void SetMonochromeMode(RDKit::MolDraw2D {lvalue},boost::python::tuple,boost::python::tuple)
    """
    ...

@overload
def SetMonochromeMode(
    self, drawer: MolDraw2D, fgColour: tuple, bgColour: tuple
) -> None: ...
def UpdateDrawerParamsFromJSON(self, drawer: MolDraw2D, json: str) -> None:
    """
    C++ signature :void UpdateDrawerParamsFromJSON(RDKit::MolDraw2D {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def UpdateMolDrawOptionsFromJSON(self, opts: MolDrawOptions, json: str) -> None:
    """
    C++ signature :void UpdateMolDrawOptionsFromJSON(RDKit::MolDrawOptions {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...
