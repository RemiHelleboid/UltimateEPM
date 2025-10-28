# EPM PARAMETERS FILE

from dataclasses import dataclass

@dataclass
class MaterialParameterEPM:
    name: str
    lattice_constant: float
    symbol: str
    V3S: float
    V3S: float
    V8S: float
    V11S: float
    V3A: float
    V4A: float
    V11A: float

# Define the parameters for the materials (Silicon, etc.)

Si = MaterialParameterEPM(
    name="Silicon",
    symbol="Si",
    lattice_constant=5.43,
    V3S=-0.263,
    V8S=-0.040,
    V11S=0.033,
    V3A=0.0,
    V4A=0.0,
    V11A=0.0,
)

# Si = MaterialParameterEPM(
#     name="Silicon",
#     symbol="Si",
#     lattice_constant=5.43,
#     V3S=-0.2241,
#     V8S=0.0551,
#     V11S=0.0724,
#     V3A=0.0,
#     V4A=0.0,
#     V11A=0.0,
# )

# Si = MaterialParameterEPM(
#     name="Silicon",
#     symbol="Si",
#     lattice_constant=5.43,
#     V3S=-0.2241,
#     V8S=0.0551,
#     V11S=0.0724,
#     V3A=0.0,
#     V4A=0.0,
#     V11A=0.0,
# )

Ge = MaterialParameterEPM(
    name="Germanium",
    symbol="Ge",
    lattice_constant=5.66,
    V3S=-0.23,
    V8S=0.01,
    V11S=0.06,
    V3A=0.0,
    V4A=0.0,
    V11A=0.0,
)

Sn = MaterialParameterEPM(
    name="Tin",
    symbol="Sn",
    lattice_constant=6.49,
    V3S=-0.20,
    V8S=0.00,
    V11S=0.04,
    V3A=0.0,
    V4A=0.0,
    V11A=0.0,
)

GaP = MaterialParameterEPM(
    name="GalliumPhosphide",
    symbol="GaP",
    lattice_constant=5.44,
    V3S=-0.22,
    V8S=0.03,
    V11S=0.07,
    V3A=0.12,
    V4A=0.07,
    V11A=0.02,
)

GaAs = MaterialParameterEPM(
    name="GalliumArsenide",
    symbol="GaAs",
    lattice_constant=5.64,
    V3S=-0.23,
    V8S=0.01,
    V11S=0.06,
    V3A=0.07,
    V4A=0.05,
    V11A=0.01,
)

AlSb = MaterialParameterEPM(
    name="AluminiumAntimonide",
    symbol="AlSb",
    lattice_constant=6.13,
    V3S=-0.21,
    V8S=0.02,
    V11S=0.06,
    V3A=0.06,
    V4A=0.04,
    V11A=0.02,
)

InP = MaterialParameterEPM(
    name="IndiumPhosphide",
    symbol="InP",
    lattice_constant=5.86,
    V3S=-0.23,
    V8S=0.01,
    V11S=0.06,
    V3A=0.07,
    V4A=0.05,
    V11A=0
)

InAs = MaterialParameterEPM(
    name="IndiumArsenide",
    symbol="InAs",
    lattice_constant=6.12,
    V3S=-0.22,
    V8S=0.00,
    V11S=0.05,
    V3A=0.08,
    V4A=0.05,
    V11A=0.03,
)

InSb = MaterialParameterEPM(
    name="IndiumAntimonide",
    symbol="InSb",
    lattice_constant=6.48,
    V3S=-0.20,
    V8S=0.00,
    V11S=0.04,
    V3A=0.06,
    V4A=0.05,
    V11A=0.01,
)

AlAs = MaterialParameterEPM(
    name="AluminumArsenide",
    symbol="AlAs",
    lattice_constant=5.6605,
    V3S=-0.221,
    V8S=0.025,
    V11S=0.07,
    V3A=0.08,
    V4A=0.05,
    V11A=-0.04,
)

ALL_MATERIALS = [Si, Ge, Sn, GaP, GaAs, AlSb, InP, InAs, InSb, AlAs]

def get_material_by_symbol(name_or_symbol):
    for material in ALL_MATERIALS:
        if material.symbol == name_or_symbol or material.name == name_or_symbol:
            return material
    raise ValueError(f"Material with name or symbol '{name_or_symbol}' not found")
