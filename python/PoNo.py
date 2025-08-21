import gmsh
import sys

gmsh.initialize()
gmsh.model.add("IBZ_Wedge")

# -----------------------------------------------------------------
# Mesh size
# -----------------------------------------------------------------
h = 0.01
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", h)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", h)

# -----------------------------------------------------------------
# Corner points of the irreducible wedge (units 2π/a)
# -----------------------------------------------------------------
p1 = gmsh.model.occ.addPoint(0.0, 0.0, 0.0, h)   # Γ
p2 = gmsh.model.occ.addPoint(0.5, 0.5, 0.5, h)   # L
p3 = gmsh.model.occ.addPoint(1.0, 0.0, 0.0, h)   # X
p4 = gmsh.model.occ.addPoint(0.75, 0.75, 0.0, h) # K
p5 = gmsh.model.occ.addPoint(1.0, 0.5, 0.0, h)   # W
p6 = gmsh.model.occ.addPoint(1.0, 0.25, 0.25, h) # U

# -----------------------------------------------------------------
# KLUW face (quad)
# -----------------------------------------------------------------
l1 = gmsh.model.occ.addLine(p4, p2)  # K-L
l2 = gmsh.model.occ.addLine(p2, p6)  # L-U
l3 = gmsh.model.occ.addLine(p6, p5)  # U-W
l4 = gmsh.model.occ.addLine(p5, p4)  # W-K
ll10 = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
s10 = gmsh.model.occ.addPlaneSurface([ll10])

# -----------------------------------------------------------------
# UWX face (triangle)
# -----------------------------------------------------------------
l6 = gmsh.model.occ.addLine(p5, p3)  # W-X
l7 = gmsh.model.occ.addLine(p3, p6)  # X-U
ll11 = gmsh.model.occ.addCurveLoop([l3, l6, l7])
s11 = gmsh.model.occ.addPlaneSurface([ll11])

# -----------------------------------------------------------------
# ΓLK face (triangle)
# -----------------------------------------------------------------
l8  = gmsh.model.occ.addLine(p1, p2) # Γ-L
l10 = gmsh.model.occ.addLine(p4, p1) # K-Γ
ll12 = gmsh.model.occ.addCurveLoop([l8, -l1, l10])
s12 = gmsh.model.occ.addPlaneSurface([ll12])

# -----------------------------------------------------------------
# ΓLUX face (quad) (Γ-L-U-X)
# -----------------------------------------------------------------
l12 = gmsh.model.occ.addLine(p3, p1) # X-Γ
ll13 = gmsh.model.occ.addCurveLoop([l8, l2, -l7, l12])
s13 = gmsh.model.occ.addPlaneSurface([ll13])

# -----------------------------------------------------------------
# ΓKWX face (quad) (Γ-K-W-X)
# -----------------------------------------------------------------
ll14 = gmsh.model.occ.addCurveLoop([-l10, -l4, l6, l12])
s14 = gmsh.model.occ.addPlaneSurface([ll14])

# -----------------------------------------------------------------
# Close the IBZ wedge into a volume
# -----------------------------------------------------------------
sl20 = gmsh.model.occ.addSurfaceLoop([s10, s11, s12, s13, s14])
v21  = gmsh.model.occ.addVolume([sl20])

gmsh.model.occ.synchronize()

# -----------------------------------------------------------------
# Physical groups for export
# -----------------------------------------------------------------
gmsh.model.addPhysicalGroup(3, [v21], 1)
gmsh.model.setPhysicalName(3, 1, "IBZ_Wedge")

gmsh.model.addPhysicalGroup(2, [s10, s11, s12, s13, s14], 2)
gmsh.model.setPhysicalName(2, 2, "IBZ_Faces")

# -----------------------------------------------------------------
# Mesh
# -----------------------------------------------------------------
gmsh.model.mesh.generate(3)
gmsh.write("ibz.msh")

if "-nopopup" not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
