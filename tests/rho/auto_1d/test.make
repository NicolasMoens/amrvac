SETUP_FLAGS := -d=1
SCHEME_DIR := ../../schemes

TESTS := ball_1d_1step_tvd.log ball_1d_2step_tvdlf_mm.log			\
ball_1d_2step_tvdmu_al.log ball_1d_3step_hll_cada.log ball_1d_4step_hll_mc.log	\
ball_1d_4step_hllc_ko.log ball_1d_rk4_tvdlf_cada.log ball_1d_ssprk43_fd_mp5.log	\
ball_1d_ssprk54_fd_mp5.log ball_1d_trap_hll_ko.log ball_1d_rev_1step_tvd.log	\
ball_1d_rev_2step_tvdlf_mm.log ball_1d_rev_2step_tvdmu_al.log			\
ball_1d_rev_3step_hll_cada.log ball_1d_rev_4step_hll_mc.log			\
ball_1d_rev_4step_hllc_ko.log ball_1d_rev_rk4_tvdlf_cada.log			\
ball_1d_rev_ssprk43_fd_mp5.log ball_1d_rev_ssprk54_fd_mp5.log			\
ball_1d_rev_trap_hll_ko.log

PARS := ball_1d.par

ball_1d_1step_tvd.log: PARS += $(SCHEME_DIR)/1step_tvd.par
ball_1d_2step_tvdlf_mm.log: PARS += $(SCHEME_DIR)/2step_tvdlf_mm.par
ball_1d_2step_tvdmu_al.log: PARS += $(SCHEME_DIR)/2step_tvdmu_al.par
ball_1d_trap_hll_ko.log: PARS += $(SCHEME_DIR)/trap_hll_ko.par
ball_1d_3step_hll_cada.log: PARS += $(SCHEME_DIR)/3step_hll_cada.par
ball_1d_4step_hll_mc.log: PARS += $(SCHEME_DIR)/4step_hll_mc.par
ball_1d_4step_hllc_ko.log: PARS += $(SCHEME_DIR)/4step_hllc_ko.par
ball_1d_rk4_tvdlf_cada.log: PARS += $(SCHEME_DIR)/rk4_tvdlf_cada.par
ball_1d_ssprk43_fd_mp5.log: PARS += $(SCHEME_DIR)/ssprk43_fd_mp5.par
ball_1d_ssprk54_fd_mp5.log: PARS += $(SCHEME_DIR)/ssprk54_fd_mp5.par

ball_1d_rev_1step_tvd.log: PARS += rev.par $(SCHEME_DIR)/1step_tvd.par
ball_1d_rev_2step_tvdlf_mm.log: PARS += rev.par $(SCHEME_DIR)/2step_tvdlf_mm.par
ball_1d_rev_2step_tvdmu_al.log: PARS += rev.par $(SCHEME_DIR)/2step_tvdmu_al.par
ball_1d_rev_trap_hll_ko.log: PARS += rev.par $(SCHEME_DIR)/trap_hll_ko.par
ball_1d_rev_3step_hll_cada.log: PARS += rev.par $(SCHEME_DIR)/3step_hll_cada.par
ball_1d_rev_4step_hll_mc.log: PARS += rev.par $(SCHEME_DIR)/4step_hll_mc.par
ball_1d_rev_4step_hllc_ko.log: PARS += rev.par $(SCHEME_DIR)/4step_hllc_ko.par
ball_1d_rev_rk4_tvdlf_cada.log: PARS += rev.par $(SCHEME_DIR)/rk4_tvdlf_cada.par
ball_1d_rev_ssprk43_fd_mp5.log: PARS += rev.par $(SCHEME_DIR)/ssprk43_fd_mp5.par
ball_1d_rev_ssprk54_fd_mp5.log: PARS += rev.par $(SCHEME_DIR)/ssprk54_fd_mp5.par

include ../../test_rules.make
