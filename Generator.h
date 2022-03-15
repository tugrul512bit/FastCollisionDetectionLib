/*
 * Generator.h
 *
 *  Created on: Mar 8, 2022
 *      Author: tugrul
 */

#ifndef GENERATOR_H_
#define GENERATOR_H_

#include<memory>

namespace oofrng
{



	template<int LANES=64>
	class Generator
	{
	public:
		Generator()
		{
			alignedSeedBuf = std::make_shared<AlignedSeedBuffer>();
		}


		const uint32_t generate1()
		{
			return rnd(alignedSeedBuf->seed);
		}

		const uint32_t generate1(const uint32_t limit)
		{
			return rnd(alignedSeedBuf->seed,limit);
		}

		const float generate1Float()
		{
			return ((float)rnd(alignedSeedBuf->seed))*alignedSeedBuf->multiplier;
		}

		const float generate1Float(const float limit)
		{
			return ((float)rnd(alignedSeedBuf->seed))*alignedSeedBuf->multiplier*limit;
		}

		// fills array of length n with values between 0 and max(2^32-1)
		void generate(uint32_t * const __restrict__ out, const size_t n)
		{

			const size_t nL = n-n%LANES;
			for(size_t i=0;i<nL;i+=LANES)
			{
				rndL(alignedSeedBuf->ptrL,out+i);
			}

			for(size_t i=nL;i<n;i++)
			{
				out[i]=rnd(alignedSeedBuf->seed);
			}
		}

		// fills array of length n with values in range [0,limit)
		void generate(uint32_t * const __restrict__ out, const size_t n, const uint32_t limit)
		{

			const size_t nL = n-n%LANES;
			for(size_t i=0;i<nL;i+=LANES)
			{
				rndL(alignedSeedBuf->ptrL,out+i,limit);
			}

			for(size_t i=nL;i<n;i++)
			{
				out[i]=rnd(alignedSeedBuf->seed,limit);
			}
		}

		// generate [0,1)
		void generate(float * const __restrict__ out, const size_t n)
		{

			const size_t nL = n-n%LANES;
			for(size_t i=0;i<nL;i+=LANES)
			{
				rndL(alignedSeedBuf->ptrL,out+i);
			}

			for(size_t i=nL;i<n;i++)
			{
				out[i]=rnd(alignedSeedBuf->seed);
			}
		}

		// generate [0,limit)
		void generate(float * const __restrict__ out, const size_t n, const float limit)
		{

			const size_t nL = n-n%LANES;
			for(size_t i=0;i<nL;i+=LANES)
			{
				rndL(alignedSeedBuf->ptrL,out+i,limit);
			}

			for(size_t i=nL;i<n;i++)
			{
				out[i]=rnd(alignedSeedBuf->seed,limit);
			}
		}

	private:

		static uint32_t* computeAlignment(uint32_t* ptr)
		{
			uint32_t* ptrLTmp = ptr;

			while(((size_t)ptrLTmp)%4096 != 0)
			{
				ptrLTmp++;
			}
			return ptrLTmp;
		}


		class AlignedSeedBuffer
		{
		public:
			AlignedSeedBuffer():ptrL(computeAlignment(seedL)),uint32_tmax(((uint32_t)0)-1),multiplier(1.0/uint32_tmax)
			{
				for(size_t i=0;i<LANES;i++)
				{
					ptrL[i]=i+1;
				}
				seed=LANES+1;
			}


			uint32_t seedL[LANES+4096];
			uint32_t seed;
			uint32_t* const __restrict__ ptrL;
			const uint32_t uint32_tmax;
			const float multiplier;
		};

		std::shared_ptr<AlignedSeedBuffer> alignedSeedBuf;



		// generate random number in range [0,max)
		const uint32_t rnd(uint32_t& seed)
		{
			// Thomas Wang's invention
			seed = (seed ^ 61) ^ (seed >> 16);
			seed *= 9;
			seed = seed ^ (seed >> 4);
			seed *= 0x27d4eb2d;
			seed = seed ^ (seed >> 15);
			return seed;
		}

		// generate random number in range [0,limit)
		const uint32_t rnd(uint32_t& seed, const uint32_t limit)
		{
			// Thomas Wang's invention
			seed = (seed ^ 61) ^ (seed >> 16);
			seed *= 9;
			seed = seed ^ (seed >> 4);
			seed *= 0x27d4eb2d;
			seed = seed ^ (seed >> 15);
			return seed%limit;
		}

		// generate [0,max)
		inline
		void rndL(uint32_t * const __restrict__ seed, uint32_t * const __restrict__ out)
		{


			for(int i=0;i<LANES;i+=2)
			{
			   const uint32_t sd = seed[i];
			   const uint32_t sd_ = seed[i+1];

			   const uint32_t sd2 = (sd ^ 61) ^ (sd >> 16);
			   const uint32_t sd2_ = (sd_ ^ 61) ^ (sd_ >> 16);

			   const uint32_t sd3 = sd2*9;
			   const uint32_t sd3_ = sd2_*9;

			   const uint32_t sd4 = sd3 ^ (sd3 >> 4);
			   const uint32_t sd4_ = sd3_ ^ (sd3_ >> 4);

			   const uint32_t sd5 = sd4*0x27d4eb2d;
			   const uint32_t sd5_ = sd4_*0x27d4eb2d;

			   const uint32_t sd6 = sd5 ^ (sd5 >> 15);
			   const uint32_t sd6_ = sd5_ ^ (sd5_ >> 15);


			   out[i]=sd6;
			   out[i+1]=sd6_;

			   seed[i]=sd6;
			   seed[i+1]=sd6_;

			}
		}

		// generate [0,limit)
		inline
		void rndL(uint32_t * const __restrict__ seed, uint32_t * const __restrict__ out, const uint32_t limit)
		{


			for(int i=0;i<LANES;i+=2)
			{
			   const uint32_t sd = seed[i];
			   const uint32_t sd_ = seed[i+1];

			   const uint32_t sd2 = (sd ^ 61) ^ (sd >> 16);
			   const uint32_t sd2_ = (sd_ ^ 61) ^ (sd_ >> 16);

			   const uint32_t sd3 = sd2*9;
			   const uint32_t sd3_ = sd2_*9;

			   const uint32_t sd4 = sd3 ^ (sd3 >> 4);
			   const uint32_t sd4_ = sd3_ ^ (sd3_ >> 4);

			   const uint32_t sd5 = sd4*0x27d4eb2d;
			   const uint32_t sd5_ = sd4_*0x27d4eb2d;

			   const uint32_t sd6 = sd5 ^ (sd5 >> 15);
			   const uint32_t sd6_ = sd5_ ^ (sd5_ >> 15);


			   out[i]=sd6%limit;
			   out[i+1]=sd6_%limit;

			   seed[i]=sd6;
			   seed[i+1]=sd6_;

			}
		}

		// generate [0,1)
		inline
		void rndL(uint32_t * const __restrict__ seed, float * const __restrict__ out)
		{

			const float mult = alignedSeedBuf->multiplier;
			for(int i=0;i<LANES;i+=2)
			{


			   const uint32_t sd = seed[i];
			   const uint32_t sd_ = seed[i+1];

			   const uint32_t sd2 = (sd ^ 61) ^ (sd >> 16);
			   const uint32_t sd2_ = (sd_ ^ 61) ^ (sd_ >> 16);

			   const uint32_t sd3 = sd2*9;
			   const uint32_t sd3_ = sd2_*9;

			   const uint32_t sd4 = sd3 ^ (sd3 >> 4);
			   const uint32_t sd4_ = sd3_ ^ (sd3_ >> 4);

			   const uint32_t sd5 = sd4*0x27d4eb2d;
			   const uint32_t sd5_ = sd4_*0x27d4eb2d;

			   const uint32_t sd6 = sd5 ^ (sd5 >> 15);
			   const uint32_t sd6_ = sd5_ ^ (sd5_ >> 15);

			   const float sd7 = sd6*mult;
			   const float sd7_ = sd6_*mult;

			   out[i]=sd7;
			   out[i+1]=sd7_;

			   seed[i]=sd6;
			   seed[i+1]=sd6_;

			}
		}

		// generate [0,limit)
		inline
		void rndL(uint32_t * const __restrict__ seed, float * const __restrict__ out, const float limit)
		{

			const float mult = alignedSeedBuf->multiplier*limit;
			for(int i=0;i<LANES;i+=2)
			{


			   const uint32_t sd = seed[i];
			   const uint32_t sd_ = seed[i+1];

			   const uint32_t sd2 = (sd ^ 61) ^ (sd >> 16);
			   const uint32_t sd2_ = (sd_ ^ 61) ^ (sd_ >> 16);

			   const uint32_t sd3 = sd2*9;
			   const uint32_t sd3_ = sd2_*9;

			   const uint32_t sd4 = sd3 ^ (sd3 >> 4);
			   const uint32_t sd4_ = sd3_ ^ (sd3_ >> 4);

			   const uint32_t sd5 = sd4*0x27d4eb2d;
			   const uint32_t sd5_ = sd4_*0x27d4eb2d;

			   const uint32_t sd6 = sd5 ^ (sd5 >> 15);
			   const uint32_t sd6_ = sd5_ ^ (sd5_ >> 15);

			   const float sd7 = sd6*mult;
			   const float sd7_ = sd6_*mult;

			   out[i]=sd7;
			   out[i+1]=sd7_;

			   seed[i]=sd6;
			   seed[i+1]=sd6_;

			}
		}
	};

}

#endif /* GENERATOR_H_ */
