TEST_DIRS := rho/auto_1d rho/auto_2d rho/auto_3d hd/Riemann_1D hd/Riemann_2D	\
hd/Rayleigh_Taylor_2D hd/Rayleigh_Taylor_3D hd/Kelvin_Helmholtz_2D		\
hd/blast_wave_Cartesian_2D hd/blast_wave_polar_2D hd/blast_wave_polar_stretched_2D \
hd/Richtmyer_Meshkov_dust_2D mhd/Riemann_1.75D	\
mhd/rotor_2D mhd/wake_2.5D mhd/shock_cloud_2D mhd/Kelvin_Helmholtz_2D		\
mhd/Orszag_Tang_2D mhd/blast_wave_Cartesian_2D mhd/blast_wave_polar_2D		\
mhd/blast_wave_spherical_3D mhd/ripple_2.5D mhd/prominence_Rayleigh_Taylor_2.5D	\
mhd/solar_atmosphere_2.5D mhd/radiative_cooling_3D mhd/blast_wave_polar_stretched_2D \
mhd/tilt_instability_2D mhd/doubleGEM_2.5D mhd/Longcope_Strauss_2D mhd/blast_wave_cylindrical_3D \
mhd/blast_wave_spherical_stretched_3D mhd/convection_2D mhd/Kelvin_Helmholtz_3D \
mhd/ring_thermal_conduction_2D mhd/ring_thermal_conduction_3D mhd/solar_flare_reconnection_2.5D

# Directories with altered names (useful for cleaning)
CLEAN_TEST := $(TEST_DIRS:%=clean-%)

.PHONY: all clean $(TEST_DIRS) $(CLEAN_TEST)

all: $(TEST_DIRS)

$(TEST_DIRS):
	@$(MAKE) -C $@ -f test.make

clean:	$(CLEAN_TEST)

$(CLEAN_TEST):
	$(MAKE) -C $(@:clean-%=%) clean
