// feat/feature-mfcc.cc

// Copyright 2009-2011  Karel Vesely;  Petr Motlicek
//                2016  Johns Hopkins University (author: Daniel Povey)

// See ../../COPYING for clarification regarding multiple authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
// THIS CODE IS PROVIDED *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED
// WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE,
// MERCHANTABLITY OR NON-INFRINGEMENT.
// See the Apache 2 License for the specific language governing permissions and
// limitations under the License.


#include "feat/feature-mfcc.h"



namespace kaldi {



void MfccComputer::Compute(BaseFloat signal_raw_log_energy,
                           BaseFloat vtln_warp,
                           VectorBase<BaseFloat> *signal_frame,
                           VectorBase<BaseFloat> *feature) {
  KALDI_ASSERT(signal_frame->Dim() == opts_.frame_opts.PaddedWindowSize() &&
               feature->Dim() == this->Dim());

  const MelBanks &mel_banks = *(GetMelBanks(vtln_warp));

  if (opts_.use_energy && !opts_.raw_energy)
    signal_raw_log_energy = Log(std::max<BaseFloat>(VecVec(*signal_frame, *signal_frame),
                                     std::numeric_limits<float>::epsilon()));
  
  /*VectorBase<BaseFloat> *signal_frame_dummy;
  for (int ii = 0; ii < signal_frame->Dim(); ii++)
    {
        signal_frame->Data()[ii] = 0.0;
    }
  for(int jj = 0;jj<signal_frame->Dim(); jj++){
    signal_frame_dummy->Data()[jj] = static_cast<double>(signal_frame->Data()[jj]);
  }*/
  //Using wavelib library as source code
  /*
  vector<double> dwt_output,flag;
  vector<int> length;
  vector<double> sig;
  int level = 3;
  for(int j = 0;j<signal_frame->Dim(); j++){
    sig[j] = static_cast<double>(signal_frame->Data()[j]);
  }
  dwt(sig,level,"db2",dwt_output,flag,length);
  for(int x = 0;x<dwt_output.size(); x++){
    signal_frame->Data()[x] = static_cast<BaseFloat>(dwt_output[x]);
  }*/
 /*
  int n = signal_frame->Dim();
  double *sig = new double[n];
  for(int loop = 0;loop<512; loop++){
  //  if(loop < n){
    sig[loop] = static_cast<double>(signal_frame->Data()[loop]);
    //}
    //else{
     // sig[loop] = 0.0;
    //}
  }*/
  
if (srfft_ != NULL)  // Compute FFT using the split-radix algorithm.
    srfft_->Compute(signal_frame->Data(), true);
  else  // An alternative algorithm that works for non-powers-of-two.
    RealFft(signal_frame, true);
/*
   //db2
  double c[4] = {
     0.4829629131445341, 
     0.8365163037378079, 
     0.2241438680420133, 
    -0.1294095225512603 };
  // db8
  int i;
  int j;
  int j0;
  int j1;
  int j2;
  int j3;
  //int p = 15;
  //int q;
  //int k;
  int m;
  double *y = new double[n];
  for ( int i_blah = 0; i_blah < n; i_blah++ )
    {
        y[i_blah] = sig[i_blah];
    }
  double *z;
  z = new double[n];
  for (i = 0; i < n; i++)
    {
        z[i] = 0.0;
    }
  m = n;
  //q = (p-1)/2;
     while ( 4 <= m )
        {
        i = 0;
        for ( j = 0; j < m - 1; j = j + 2 )
            {
                //i4_wrap for j0 - sending j,0,m-1
                int jhi_0;
                int jlo_0;
                int value_0;
                int wide_0;
                int ival_0 = j;
                int ilo_0 = 0;
                int ihi_0 = m-1;
                //start i4min for j0 - sending it ilo,ihi ==jlo
                int i1_0 = ilo_0;
                int i2_0 = ihi_0;
                if ( i1_0 < i2_0 )
                {
                    jlo_0 = i1_0;
                }
                else
                {
                    jlo_0 = i2_0;
                }
                //end  i4min for j0
                //start i4max for j0 - sending it ilo,ihi == jhi
                
                if ( i2_0 < i1_0 )
                {
                    jhi_0 = i1_0;
                }
                else
                {
                    jhi_0 = i2_0;
                }
                //end i4max for j0
                wide_0 = jhi_0 + 1 - jlo_0;
                if ( wide_0 == 1 )
                    {
                        value_0 = jlo_0;
                    }
                else
                {
                    int temp_0;
                    //start of i4mop sending it ival_0-jlo_0,wide_0 == temp_0
                    int i_0 = ival_0-jlo_0;
                    int j_0 = wide_0;
                    temp_0 = i_0 % j_0;

                    if ( temp_0 < 0 )
                    {
                        temp_0 = temp_0 + abs ( j_0 );
                    }
                    //end i4modp_0
                    value_0 = jlo_0 + temp_0;
                }
                j0 = value_0;
                //end i4wrap_0
                //i4_wrap for j1 - sending j+1,0,n-1
                int jhi_1;
                int jlo_1;
                int value_1;
                int wide_1;
                int ival_1 = j+1;
                int ilo_1 = 0;
                int ihi_1 = m-1;
                //start i4min for j1 - sending it ilo,ihi ==jlo
                int i1_1 = ilo_1;
                int i2_1 = ihi_1;
                if ( i1_1 < i2_1 )
                {
                    jlo_1 = i1_1;
                }
                else
                {
                    jlo_1 = i2_1;
                }
                //end  i4min for j1
                //start i4max for j1 - sending it ilo,ihi == jhi
                
                if ( i2_1 < i1_1 )
                {
                    jhi_1 = i1_1;
                }
                else
                {
                    jhi_1 = i2_1;
                }
                //end i4max for j1
                wide_1 = jhi_1 + 1 - jlo_1;
                if ( wide_1 == 1 )
                    {
                        value_1 = jlo_1;
                    }
                else
                {
                    int temp_1;
                    //start of i4mop sending it ival_1-jlo_1,wide_0 == temp_1
                    int i_1 = ival_1-jlo_1;
                    int j_1 = wide_1;
                    temp_1 = i_1 % j_1;

                    if ( temp_1 < 0 )
                    {
                        temp_1 = temp_1 + abs ( j_1 );
                    }
                    //end i4modp_1
                    value_1 = jlo_1 + temp_1;
                }
                j1 = value_1;
                //end i4wrap_1
                
                //i4_wrap for j2 - sending j+2,0,n-1
                int jhi_2;
                int jlo_2;
                int value_2;
                int wide_2;
                int ival_2 = j+2;
                int ilo_2 = 0;
                int ihi_2 = m-1;
                //start i4min for j2 - sending it ilo,ihi ==jlo
                int i1_2 = ilo_2;
                int i2_2 = ihi_2;
                if ( i1_2 < i2_2 )
                {
                    jlo_2 = i1_2;
                }
                else
                {
                    jlo_2 = i2_2;
                }
                //end  i4min for j2
                //start i4max for j2 - sending it ilo,ihi == jhi
                
                if ( i2_2 < i1_2 )
                {
                    jhi_2 = i1_2;
                }
                else
                {
                    jhi_2 = i2_2;
                }
                //end i4max for j1
                wide_2 = jhi_2 + 1 - jlo_2;
                if ( wide_2 == 1 )
                    {
                        value_2 = jlo_2;
                    }
                else
                {
                    int temp_2;
                    //start of i4mop sending it ival_2-jlo_2,wide_0 == temp_2
                    int i_2 = ival_2-jlo_2;
                    int j_2 = wide_2;
                    temp_2 = i_2 % j_2;

                    if ( temp_2 < 0 )
                    {
                        temp_2 = temp_2 + abs ( j_2 );
                    }
                    //end i4modp_2
                    value_2 = jlo_2 + temp_2;
                }
                j2 = value_2;
                //end i4wrap_1
                //i4_wrap for j3 - sending j+3,0,n-1
                int jhi_3;
                int jlo_3;
                int value_3;
                int wide_3;
                int ival_3 = j+3;
                int ilo_3 = 0;
                int ihi_3 = m-1;
                //start i4min for j3 - sending it ilo,ihi ==jlo
                int i1_3 = ilo_3;
                int i2_3 = ihi_3;
                if ( i1_3 < i2_3 )
                {
                    jlo_3 = i1_3;
                }
                else
                {
                    jlo_3 = i2_3;
                }
                //end  i4min for j3
                //start i4max for j3 - sending it ilo,ihi == jhi
                
                if ( i2_3 < i1_3 )
                {
                    jhi_3 = i1_3;
                }
                else
                {
                    jhi_3 = i2_3;
                }
                //end i4max for j3
                wide_3 = jhi_3 + 1 - jlo_3;
                if ( wide_3 == 1 )
                    {
                        value_3 = jlo_3;
                    }
                else
                {
                    int temp_3;
                    //start of i4mop sending it ival_3-jlo_3,wide_0 == temp_3
                    int i_3 = ival_3-jlo_3;
                    int j_3 = wide_3;
                    temp_3 = i_3 % j_3;

                    if ( temp_3 < 0 )
                    {
                        temp_3 = temp_3 + abs ( j_3 );
                    }
                    //end i4modp_3
                    value_3 = jlo_3 + temp_3;
                }
                j3 = value_3;
                //end i4wrap_3
                
                z[i] = c[0] * y[j0] + c[1] * y[j1] + c[2] * y[j2] + c[3] * y[j3];
                z[i+m/2] = c[3] * y[j0] - c[2] * y[j1] + c[1] * y[j2] - c[0] * y[j3];
          
                i = i + 1;
            }
        for ( i = 0; i < m; i++ )
            {
            y[i] = z[i];
            }
        m = m / 2;
        }
    for(int b = 0; b<n;b++){
      signal_frame->Data()[b] = static_cast<BaseFloat>(y[b]);
    }

*/
  
  
  // Convert the FFT into a power spectrum.
  ComputePowerSpectrum(signal_frame);
  SubVector<BaseFloat> power_spectrum(*signal_frame, 0,
                                      signal_frame->Dim() / 2 + 1);

  //Convert vector double output into a SubVector for mel.banks Compute
  //mel_banks.Compute(dwt_output, &mel_energies_);
  mel_banks.Compute(power_spectrum, &mel_energies_);  
  
  // avoid log of zero (which should be prevented anyway by dithering).
  mel_energies_.ApplyFloor(std::numeric_limits<float>::epsilon());
  mel_energies_.ApplyLog();  // take the log.

  feature->SetZero();  // in case there were NaNs.
  // feature = dct_matrix_ * mel_energies [which now have log]
  feature->AddMatVec(1.0, dct_matrix_, kNoTrans, mel_energies_, 0.0);

  if (opts_.cepstral_lifter != 0.0)
    feature->MulElements(lifter_coeffs_);

  if (opts_.use_energy) {
    if (opts_.energy_floor > 0.0 && signal_raw_log_energy < log_energy_floor_)
      signal_raw_log_energy = log_energy_floor_;
    (*feature)(0) = signal_raw_log_energy;
  }

  if (opts_.htk_compat) {
    BaseFloat energy = (*feature)(0);
    for (int32 i = 0; i < opts_.num_ceps - 1; i++)
      (*feature)(i) = (*feature)(i+1);
    if (!opts_.use_energy)
      energy *= M_SQRT2;  // scale on C0 (actually removing a scale
    // we previously added that's part of one common definition of
    // the cosine transform.)
    (*feature)(opts_.num_ceps - 1)  = energy;
  }
}

MfccComputer::MfccComputer(const MfccOptions &opts):
    opts_(opts), srfft_(NULL),
    mel_energies_(opts.mel_opts.num_bins) {

  int32 num_bins = opts.mel_opts.num_bins;
  if (opts.num_ceps > num_bins)
    KALDI_ERR << "num-ceps cannot be larger than num-mel-bins."
              << " It should be smaller or equal. You provided num-ceps: "
              << opts.num_ceps << "  and num-mel-bins: "
              << num_bins;

  Matrix<BaseFloat> dct_matrix(num_bins, num_bins);
  ComputeDctMatrix(&dct_matrix);
  // Note that we include zeroth dct in either case.  If using the
  // energy we replace this with the energy.  This means a different
  // ordering of features than HTK.
  SubMatrix<BaseFloat> dct_rows(dct_matrix, 0, opts.num_ceps, 0, num_bins);
  dct_matrix_.Resize(opts.num_ceps, num_bins);
  dct_matrix_.CopyFromMat(dct_rows);  // subset of rows.
  if (opts.cepstral_lifter != 0.0) {
    lifter_coeffs_.Resize(opts.num_ceps);
    ComputeLifterCoeffs(opts.cepstral_lifter, &lifter_coeffs_);
  }
  if (opts.energy_floor > 0.0)
    log_energy_floor_ = Log(opts.energy_floor);

  int32 padded_window_size = opts.frame_opts.PaddedWindowSize();
  if ((padded_window_size & (padded_window_size-1)) == 0)  // Is a power of two...
    srfft_ = new SplitRadixRealFft<BaseFloat>(padded_window_size);

  // We'll definitely need the filterbanks info for VTLN warping factor 1.0.
  // [note: this call caches it.]
  GetMelBanks(1.0);
}

MfccComputer::MfccComputer(const MfccComputer &other):
    opts_(other.opts_), lifter_coeffs_(other.lifter_coeffs_),
    dct_matrix_(other.dct_matrix_),
    log_energy_floor_(other.log_energy_floor_),
    mel_banks_(other.mel_banks_),
    srfft_(NULL),
    mel_energies_(other.mel_energies_.Dim(), kUndefined) {
  for (std::map<BaseFloat, MelBanks*>::iterator iter = mel_banks_.begin();
       iter != mel_banks_.end(); ++iter)
    iter->second = new MelBanks(*(iter->second));
  if (other.srfft_ != NULL)
    srfft_ = new SplitRadixRealFft<BaseFloat>(*(other.srfft_));
}



MfccComputer::~MfccComputer() {
  for (std::map<BaseFloat, MelBanks*>::iterator iter = mel_banks_.begin();
      iter != mel_banks_.end();
      ++iter)
    delete iter->second;
  delete srfft_;
}

const MelBanks *MfccComputer::GetMelBanks(BaseFloat vtln_warp) {
  MelBanks *this_mel_banks = NULL;
  std::map<BaseFloat, MelBanks*>::iterator iter = mel_banks_.find(vtln_warp);
  if (iter == mel_banks_.end()) {
    this_mel_banks = new MelBanks(opts_.mel_opts,
                                  opts_.frame_opts,
                                  vtln_warp);
    mel_banks_[vtln_warp] = this_mel_banks;
  } else {
    this_mel_banks = iter->second;
  }
  return this_mel_banks;
}



}  // namespace kaldi
