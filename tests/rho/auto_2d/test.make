SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes

TESTS := vac_2d_1step_tvd.log vac_2d_2step_tvdlf_mm.log				\
vac_2d_2step_tvdmu_al.log vac_2d_3step_hll_cada.log vac_2d_4step_hll_mc.log	\
vac_2d_4step_hllc_ko.log vac_2d_rk4_tvdlf_cada.log vac_2d_ssprk43_fd_mp5.log	\
vac_2d_ssprk54_fd_mp5.log ball_2d_1step_tvd.log ball_2d_2step_tvdlf_mm.log	\
ball_2d_2step_tvdmu_al.log ball_2d_3step_hll_cada.log ball_2d_4step_hll_mc.log	\
ball_2d_4step_hllc_ko.log ball_2d_rk4_tvdlf_cada.log ball_2d_ssprk43_fd_mp5.log	\
ball_2d_ssprk54_fd_mp5.log ball_2d_rev_1step_tvd.log				\
ball_2d_rev_2step_tvdlf_mm.log ball_2d_rev_2step_tvdmu_al.log			\
ball_2d_rev_3step_hll_cada.log ball_2d_rev_4step_hll_mc.log			\
ball_2d_rev_4step_hllc_ko.log ball_2d_rev_rk4_tvdlf_cada.log			\
ball_2d_rev_ssprk43_fd_mp5.log ball_2d_rev_ssprk54_fd_mp5.log			\
ball_2d_trap_hll_ko.log ball_2d_rev_trap_hll_ko.log

vac_2d_1step_tvd.log: PARS = vac_2d.par $(SCHEME_DIR)/1step_tvd.par
vac_2d_2step_tvdlf_mm.log: PARS = vac_2d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
vac_2d_2step_tvdmu_al.log: PARS = vac_2d.par $(SCHEME_DIR)/2step_tvdmu_al.par
vac_2d_3step_hll_cada.log: PARS = vac_2d.par $(SCHEME_DIR)/3step_hll_cada.par
vac_2d_4step_hll_mc.log: PARS = vac_2d.par $(SCHEME_DIR)/4step_hll_mc.par
vac_2d_4step_hllc_ko.log: PARS = vac_2d.par $(SCHEME_DIR)/4step_hllc_ko.par
vac_2d_rk4_tvdlf_cada.log: PARS = vac_2d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par
vac_2d_ssprk43_fd_mp5.log: PARS = vac_2d.par $(SCHEME_DIR)/ssprk43_fd_mp5.par
vac_2d_ssprk54_fd_mp5.log: PARS = vac_2d.par $(SCHEME_DIR)/ssprk54_fd_mp5.par

ball_2d_1step_tvd.log: PARS = ball_2d.par $(SCHEME_DIR)/1step_tvd.par
ball_2d_2step_tvdlf_mm.log: PARS = ball_2d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
ball_2d_2step_tvdmu_al.log: PARS = ball_2d.par $(SCHEME_DIR)/2step_tvdmu_al.par
ball_2d_3step_hll_cada.log: PARS = ball_2d.par $(SCHEME_DIR)/3step_hll_cada.par
ball_2d_4step_hll_mc.log: PARS = ball_2d.par $(SCHEME_DIR)/4step_hll_mc.par
ball_2d_4step_hllc_ko.log: PARS = ball_2d.par $(SCHEME_DIR)/4step_hllc_ko.par
ball_2d_rk4_tvdlf_cada.log: PARS = ball_2d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par
ball_2d_ssprk43_fd_mp5.log: PARS = ball_2d.par $(SCHEME_DIR)/ssprk43_fd_mp5.par
ball_2d_ssprk54_fd_mp5.log: PARS = ball_2d.par $(SCHEME_DIR)/ssprk54_fd_mp5.par
ball_2d_trap_hll_ko.log: PARS = ball_2d.par $(SCHEME_DIR)/trap_hll_ko.par

ball_2d_rev_1step_tvd.log: PARS = ball_2d.par rev.par $(SCHEME_DIR)/1step_tvd.par
ball_2d_rev_2step_tvdlf_mm.log: PARS = ball_2d.par rev.par $(SCHEME_DIR)/2step_tvdlf_mm.par
ball_2d_rev_2step_tvdmu_al.log: PARS = ball_2d.par rev.par $(SCHEME_DIR)/2step_tvdmu_al.par
ball_2d_rev_3step_hll_cada.log: PARS = ball_2d.par rev.par $(SCHEME_DIR)/3step_hll_cada.par
ball_2d_rev_4step_hll_mc.log: PARS = ball_2d.par rev.par $(SCHEME_DIR)/4step_hll_mc.par
ball_2d_rev_4step_hllc_ko.log: PARS = ball_2d.par rev.par $(SCHEME_DIR)/4step_hllc_ko.par
ball_2d_rev_rk4_tvdlf_cada.log: PARS = ball_2d.par rev.par $(SCHEME_DIR)/rk4_tvdlf_cada.par
ball_2d_rev_ssprk43_fd_mp5.log: PARS = ball_2d.par rev.par $(SCHEME_DIR)/ssprk43_fd_mp5.par
ball_2d_rev_ssprk54_fd_mp5.log: PARS = ball_2d.par rev.par $(SCHEME_DIR)/ssprk54_fd_mp5.par
ball_2d_rev_trap_hll_ko.log: PARS = ball_2d.par rev.par $(SCHEME_DIR)/trap_hll_ko.par

include ../../test_rules.make
