/*
   Comp Eng 3DY4 (Computer Systems Integration Project)
   Department of Electrical and Computer Engineering
   McMaster University
   Ontario, Canada
*/

#include <limits.h>
#include "dy4.h"
#include "iofunc.h"
#include "filter.h"
#include "gtest/gtest.h"

namespace {

    class Convolution_Fixture: public ::testing::Test {
    public:
        const int N = 1024;    
        const int M = 101;     
        const int BLOCK_SIZE = 256; 
        const int DS_SMALL = 4;   
        const real lower_bound = -1;
        const real upper_bound = 1;
        const real EPSILON = 1e-4;

        std::vector<real> x, h, y_reference, y_test, y_block_test;

        Convolution_Fixture() {
            x.resize(N);
            h.resize(M);
            y_reference.resize(N + M - 1);
            y_test.resize(N + M - 1);
            y_block_test.resize(N); // Block processing typically produces N outputs
        }

        void SetUp() {
            generate_random_values(x, lower_bound, upper_bound);
            generate_random_values(h, lower_bound, upper_bound);
            convolveFIR_reference(y_reference, x, h);
        }
    };

    TEST_F(Convolution_Fixture, convolveFIR_inefficient_NEAR) {
        convolveFIR(y_test, x, h);

        ASSERT_EQ(y_reference.size(), y_test.size()) 
            << "Output vector sizes for convolveFIR_reference and convolveFIR are unequal";

        for (int i = 0; i < (int)y_reference.size(); ++i) {
            EXPECT_NEAR(y_reference[i], y_test[i], EPSILON) 
                << "Vectors differ at index " << i;
        }
    }

    TEST_F(Convolution_Fixture, blockConvolution_NEAR) {
        std::vector<real> state(h.size() - 1, 0.0);
        std::fill(y_block_test.begin(), y_block_test.end(), 0.0);

        for (int b = 0; b < N; b += BLOCK_SIZE) {
            int block_end = std::min(b + BLOCK_SIZE, N);
            std::vector<real> x_block(x.begin() + b, x.begin() + block_end);
            std::vector<real> y_block(block_end - b, 0.0);
            
            block_conv(y_block, h, x_block, state);
            std::copy(y_block.begin(), y_block.end(), y_block_test.begin() + b);
        }
        for (int i = 0; i < N; ++i) {
            EXPECT_NEAR(y_reference[i], y_block_test[i], EPSILON)
                << "Mismatch at index " << i;
        }
    }
   
	// test downsampler
	TEST_F(Convolution_Fixture, DownSampler_NEAR) {
		const int down_sample_coeff = 4;
		std::vector<real> state(h.size() - 1, 0.0);
		std::vector<real> y_downsampled;
		std::vector<real> y_convolution;
		convolveFIR(y_convolution, x, h);
		
		std::vector<real> y_reference_downsampled(x.size() / down_sample_coeff);
		for (int i = 0; i < y_reference_downsampled.size(); i++) {
			y_reference_downsampled[i] = y_convolution[i * down_sample_coeff];
		}
		
		downSampler(y_downsampled, x, h, down_sample_coeff, state);
		
		ASSERT_EQ(y_reference_downsampled.size(), y_downsampled.size());
		for (int i = 0; i < (int)y_reference_downsampled.size(); ++i) {
			EXPECT_NEAR(y_reference_downsampled[i], y_downsampled[i], EPSILON);
		}
		for (int n = 0; n < h.size() - 1; n++) {
			int x_index = (x.size() - h.size()) + n;
			EXPECT_NEAR(x[x_index], state[n], EPSILON) 
				<< "State mismatch at index " << n;
		}
	}
} // end of namespace

