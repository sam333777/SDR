

#include <limits.h>
#include "dy4.h"
#include "iofunc.h"
#include "filter.h"
#include "gtest/gtest.h"
#include <numeric>

namespace {

    class DownSampler_Fixture: public ::testing::Test {
    public:
        const int N = 1024;
        const int M = 101;
        const int decim = 10; 
        const int lower_bound = -1;
        const int upper_bound = 1;
        const real EPSILON = 1e-4;

        std::vector<real> x, h, y_reference, y_test, state;
        

        const real Fs = 48000.0;
        const real Fc = 3000.0;
        const unsigned short int num_taps = 101;
        const int up_sample_coeff = 4;
        std::vector<real> h_standard, h_resampled;
        
        DownSampler_Fixture() {
            x.resize(N);
            h.resize(M);
            y_reference.resize(N + M - 1);
            y_test.resize(N + M - 1);
            state.resize(M - 1, 0.0);

            h_standard.resize(num_taps);
            h_resampled.resize(num_taps * up_sample_coeff);
        }
        
        void SetUp() {
            generate_random_values(x, lower_bound, upper_bound);
            generate_random_values(h, lower_bound, upper_bound);
            convolveFIR_reference(y_reference, x, h);
        }
    };

    

    TEST_F(DownSampler_Fixture, ResamplerComparison_NEAR) {
        const int up_sample_coeff = 2;
        const int down_sample_coeff = 3;
        
        std::vector<real> state_ref(M - 1, 0.0);
        std::vector<real> state_opt(M - 1, 0.0);
        std::vector<real> y_ref, y_opt;
        resampler(y_ref, x, h, down_sample_coeff, up_sample_coeff, state_ref);
        resampler_optimized(y_opt, x, h, down_sample_coeff, up_sample_coeff, state_opt);
        
        
        ASSERT_EQ(y_ref.size(), y_opt.size()) 
            << "Output sizes differ between resampler implementations";
        for (int i = 0; i < (int)y_ref.size(); ++i) {
            EXPECT_NEAR(y_ref[i], y_opt[i], EPSILON) 
                << "Resampler outputs differ at index " << i;
        }
        ASSERT_EQ(state_ref.size(), state_opt.size()) 
            << "State sizes differ between implementations";
        for (int i = 0; i < (int)state_ref.size(); ++i) {
            EXPECT_NEAR(state_ref[i], state_opt[i], EPSILON)
                << "State vectors differ at index " << i;
        }
    }
    
    TEST_F(DownSampler_Fixture, LPF_NEAR) {
        const real TOLERANCE = 0.01;
        
        impulseResponseLPF(Fs, Fc, num_taps, h_standard);
        impulseResponseLPFwithResamplerFactor(Fs, Fc, num_taps, h_resampled, up_sample_coeff);
    
        std::vector<real> resampled_decimated(num_taps);
        for (int i = 0; i < num_taps; i++) {
            int upsampled_pos = i * up_sample_coeff;
            resampled_decimated[i] = h_resampled[upsampled_pos];
        }
    
        for (int i = 0; i < num_taps; i++) {
            EXPECT_NEAR(resampled_decimated[i], h_standard[i], TOLERANCE)
                << "Mismatch at tap position " << i;
        }
    }

    TEST_F(DownSampler_Fixture, AllPassFilter_NEAR) {
        const int delay_length = state.size();
        std::vector<real> y_output;
        std::vector<real> filter_state = state;
    
        allPassFilter(y_output, x, filter_state);
    
        ASSERT_EQ(y_output.size(), x.size());
        
        for (int i = 0; i < delay_length; ++i) {
            EXPECT_NEAR(y_output[i], state[i], EPSILON);
        }
        
        for (int i = 0; i < x.size() - delay_length; ++i) {
            EXPECT_NEAR(y_output[i + delay_length], x[i], EPSILON);
        }
        
        ASSERT_EQ(filter_state.size(), delay_length);
        
        for (int i = 0; i < delay_length; ++i) {
            EXPECT_NEAR(filter_state[i], x[x.size() - delay_length + i], EPSILON);
        }
        
        std::vector<real> y_output2;
        std::vector<real> x2(x.size(), 0.0);
        generate_random_values(x2, lower_bound, upper_bound);
        
        allPassFilter(y_output2, x2, filter_state);
        
        for (int i = 0; i < delay_length; ++i) {
            EXPECT_NEAR(y_output2[i], x[x.size() - delay_length + i], EPSILON);
        }
        for (int i = 0; i < delay_length; ++i) {
            EXPECT_NEAR(filter_state[i], x2[x2.size() - delay_length + i], EPSILON);
        }
    }
    
    


} // end of namespace

