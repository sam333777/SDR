#include <limits.h>
#include <iostream>
#include "dy4.h"
#include "iofunc.h"
#include "filter.h"
#include "gtest/gtest.h"

namespace {

class PLLStateSaving_Fixture : public ::testing::Test {
public:
    const int num_samples = 1024;
    const int block_size = 128;
    const real epsilon = 1e-4;
    const real freq_tolerance = 0.1; 
    
    const real sampling_rate = 240e3;
    const real target_freq = 19e3;
    const real nco_scale = 2.0;
    const real phase_adjust = 0.0;
    const real norm_bw = 0.01;

    const int lower_bound = -1;
    const int upper_bound = 1;

    std::vector<real> x;
    std::vector<real> pll_ref, pll_test, pll_testRDS;

    PLLStateSaving_Fixture() {
        x.resize(num_samples);
        pll_ref.resize(num_samples+1);
    }

    void SetUp() override {
        generate_random_values(x, lower_bound, upper_bound);
        
        real nextncoOut = 1.0;
        fmPll(x, target_freq, sampling_rate, nextncoOut, pll_ref, 
              nco_scale, phase_adjust, norm_bw);
    }
};

TEST_F(PLLStateSaving_Fixture, PLL_block_NEAR) { 
    std::vector<real> xb, pll_b;
    real saveI = 1.0;
    real saveQ = 0.0;
    real integrator = 0.0;
    real eastmatePhase = 0.0;
    real trigOffset = 0;
    real ncoFirst = 1.0;
    
    for (int i = 0; i < num_samples; i += block_size) {
        xb.clear();
        int current_block_size = std::min(block_size, num_samples - i);
        xb.assign(x.begin() + i, x.begin() + i + current_block_size);
        pll_b.resize(current_block_size);
        
        fmPllwithState(xb, target_freq, sampling_rate, nco_scale, phase_adjust, norm_bw, saveI, saveQ, integrator, eastmatePhase, trigOffset, pll_b, ncoFirst);
        
        pll_test.insert(pll_test.end(), pll_b.begin(), pll_b.end());
    }
    ASSERT_EQ(num_samples, pll_test.size()) << " size mismatch";
    ASSERT_EQ(num_samples, pll_ref.size() - 1) << "size mismatch";
    
    for (int i = 0; i < num_samples; ++i) {
        EXPECT_NEAR(pll_ref[i], pll_test[i], epsilon)
            << "mismatch at " << i;
    }
}

}
