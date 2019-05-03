functions {
    real intensity_log_std(real z, real scaleHi, real scaleLo, real offset, real bend, real smooth) {
        return 0.5*(scaleHi+scaleLo)*(z-bend) + 0.5*(scaleHi-scaleLo)*sqrt((z-bend)*(z-bend)+smooth) + offset;
    }
}

data {
  int<lower=1> Nobservations;   // number of MS observations of pepmodXcharge
  int<lower=1> Nmschannels;     // number of MS channels (MS run X sample)
  int<lower=1> Nmsruns;         // number of MS runs
  int<lower=1> Nsamples;        // number of biological samples (if not 1-to-1 with MS runs due to labeling)
  int<lower=1> Npepmods;        // number of peptides with modifications
  int<lower=1> NpepmodCharges;  // number of peptides with modifications and charge states

  int<lower=0, upper=1> useSamples;
  int<lower=0, upper=1> useSampleEffects;
  int<lower=0, upper=1> useMsrunEffects;
  int<lower=0, upper=1> useCharges;       // if pepmod charges should be used

  int<lower=1,upper=Nmschannels*NpepmodCharges> obs2mschanXpmc[Nobservations];

  int<lower=1,upper=NpepmodCharges> mschanXpmc2pmc[Nmschannels*NpepmodCharges];
  int<lower=1,upper=Nmschannels> mschanXpmc2mschan[Nmschannels*NpepmodCharges];
  int<lower=1,upper=Npepmods> pmc2pepmod[NpepmodCharges];

  int<lower=1,upper=Nmsruns> mschan2msrun[Nmschannels];
  int<lower=1,upper=Nsamples> mschan2sample[Nmschannels];

  int<lower=0> NmsrunEffects;
  int<lower=0> NsampleEffects;
  int<lower=0> NpepmodEffects;
  int<lower=0> NpmcEffects;
  int<lower=0> NmsrunXpepmodEffects;
  int<lower=0> NsampleXpepmodEffects;

  matrix[Nmsruns, NmsrunEffects] msrun_eff2shift;
  matrix[Nsamples, NsampleEffects] sample_eff2shift;
  matrix[Npepmods, NpepmodEffects] pepmod_eff2shift;
  matrix[NpepmodCharges, NpmcEffects] pmc_eff2shift;
  matrix[Nmsruns*Npepmods, NmsrunXpepmodEffects] msrunXpepmod_eff2shift;
  matrix[Nsamples*Npepmods, NsampleXpepmodEffects] sampleXpepmod_eff2shift;

  real<lower=0.0> msrun_shift_scale_tau;
  real<lower=0.0> sample_shift_scale_tau;
  real<lower=0.0> pepmod_shift_scale_tau;
  real<lower=0.0> pmc_shift_scale_tau;
  real<lower=0.0> msrunXpepmod_scale_tau;
  real<lower=0.0> sampleXpepmod_scale_tau;

  real<lower=1.0> msrun_shift_scale_chi;
  real<lower=1.0> sample_shift_scale_chi;
  real<lower=1.0> pepmod_shift_scale_chi;
  real<lower=1.0> pmc_shift_scale_chi;
  real<lower=1.0> msrunXpepmod_scale_chi;
  real<lower=1.0> sampleXpepmod_scale_chi;

  // map from labelXreplicateXobject to observed/missed data
  int<lower=0> Nquanted;        // total number of quantified objectsXexperiments
  int<lower=1,upper=Nobservations>  quant2obs[Nquanted];
  int<lower=0> Nmissed;         // total number of missed objectsXexperiments
  int<lower=1,upper=Nobservations> miss2obs[Nmissed];

  vector<lower=0>[Nquanted] qData; // quanted data

  // global model constants
  real global_labu_shift;   // shift to be applied to all XXX_labu variables to get the real log intensity
  vector[Nmsruns] global_msrun_shifts;
  vector[Nsamples] global_sample_shifts;

  // instrument calibrated parameters 
  real<lower=0> zDetectionFactor;
  real zDetectionIntercept;
  real<lower=0, upper=1> detectionMax;

  real<lower=0> sigmaScaleHi;
  real<lower=0> sigmaScaleLo;
  real sigmaOffset;
  real sigmaBend;
  real sigmaSmooth;

  real zShift;
  real zScale;
}

transformed data {
    vector[Nquanted] qShifts;
    vector[Nmissed] mShifts;

    real mzShift; // zShift for the missing observation intensity (zShift shifted by obj_base)
    vector[Nquanted] zScore; // log(qData) transformed in zScore
    vector[Nquanted] qLogStd; // log(sd(qData))-obj_base
    vector<lower=0>[Nquanted] qDataNorm; // qData/sd(qData)

    int<lower=1,upper=Nmsruns> obs2msrun[Nobservations];
    int<lower=1,upper=Nsamples> obs2sample[Nobservations];
    int<lower=1,upper=Npepmods> obs2pepmod[Nobservations];
    int<lower=1,upper=NpepmodCharges> obs2pmc[Nobservations];
    int<lower=1,upper=Nmsruns*Npepmods> obs2msrunXpepmod[Nobservations];
    int<lower=1,upper=Nsamples*Npepmods> obs2sampleXpepmod[Nobservations];

    obs2msrun = mschan2msrun[mschanXpmc2mschan[obs2mschanXpmc]];
    obs2sample = mschan2sample[mschanXpmc2mschan[obs2mschanXpmc]];
    obs2pmc = mschanXpmc2pmc[obs2mschanXpmc];
    obs2pepmod = pmc2pepmod[obs2pmc];
    for (i in 1:Nobservations) {
      obs2msrunXpepmod[i] = (obs2msrun[i] - 1).*Npepmods + obs2pepmod[i];
      obs2sampleXpepmod[i] = (obs2sample[i] - 1).*Npepmods + obs2pepmod[i];
    }

    {
        vector[Nquanted] qLogData;
        qLogData = log(qData);
        zScore = (qLogData - zShift) * zScale;
        mzShift = zShift - global_labu_shift;
  
        // process the intensity data to optimize likelihood calculation
        for (i in 1:Nquanted) {
            qLogStd[i] = intensity_log_std(zScore[i], sigmaScaleHi, sigmaScaleLo, sigmaOffset, sigmaBend, sigmaSmooth);
            qDataNorm[i] = exp(qLogData[i] - qLogStd[i]);
            qLogStd[i] = qLogStd[i] - global_labu_shift; // obs_labu is modeled without obj_base
        }
    }

    qShifts = rep_vector(0.0, Nquanted);
    mShifts = rep_vector(0.0, Nmissed);
    if (Nmsruns > 0) {
        qShifts = qShifts + global_msrun_shifts[obs2msrun[quant2obs]];
        mShifts = mShifts + global_msrun_shifts[obs2msrun[miss2obs]];
    }
    if (Nsamples > 0) {
        qShifts = qShifts + global_sample_shifts[obs2sample[quant2obs]];
        mShifts = mShifts + global_sample_shifts[obs2sample[miss2obs]];
    }
}

parameters {
    real base_labu;

    vector[NmsrunEffects] msrun_effect;
    vector[NsampleEffects] sample_effect;
    vector[NpepmodEffects] pepmod_effect;
    vector[NpmcEffects] pmc_effect;
    vector[NmsrunXpepmodEffects] msrunXpepmod_effect;
    vector[NsampleXpepmodEffects] sampleXpepmod_effect;

    real<lower=0.0> msrun_shift_scale_lambda_t;
    real<lower=0.0> msrun_shift_scale_lambda_a;
    real<lower=0.0> sample_shift_scale_lambda_t;
    real<lower=0.0> sample_shift_scale_lambda_a;
    real<lower=0.0> pepmod_shift_scale_lambda_t;
    real<lower=0.0> pepmod_shift_scale_lambda_a;
    real<lower=0.0> pmc_shift_scale_lambda_t;
    real<lower=0.0> pmc_shift_scale_lambda_a;

    vector<lower=0.0>[NmsrunXpepmodEffects > 0 ? Nmsruns*Npepmods : 0] msrunXpepmod_shift_scale_lambda_t;
    vector<lower=0.0>[NmsrunXpepmodEffects > 0 ? Nmsruns*Npepmods : 0] msrunXpepmod_shift_scale_lambda_a;
    vector<lower=0.0>[NsampleXpepmodEffects > 0 ? Nsamples*Npepmods : 0] sampleXpepmod_shift_scale_lambda_t;
    vector<lower=0.0>[NsampleXpepmodEffects > 0 ? Nsamples*Npepmods : 0] sampleXpepmod_shift_scale_lambda_a;
}

transformed parameters {
    vector[Nobservations] o_labu;

    real<lower=0.0> msrun_shift_scale;
    real<lower=0.0> sample_shift_scale;
    real<lower=0.0> pepmod_shift_scale;
    real<lower=0.0> pmc_shift_scale;
    vector<lower=0.0>[NmsrunXpepmodEffects > 0 ? Nmsruns*Npepmods : 0] msrunXpepmod_shift_scale;
    vector<lower=0.0>[NsampleXpepmodEffects > 0 ? Nsamples*Npepmods : 0] sampleXpepmod_shift_scale;

    vector[Nmsruns] msrun_shift_unscaled;
    vector[Nsamples] sample_shift_unscaled;
    vector[Npepmods] pepmod_shift_unscaled;
    vector[NpepmodCharges] pmc_shift_unscaled;
    vector[Nmsruns*Npepmods] msrunXpepmod_shift_unscaled;
    vector[Nmsruns*Npepmods] sampleXpepmod_shift_unscaled;

    vector[Nmsruns] msrun_shift;
    vector[Nsamples] sample_shift;
    vector[Npepmods] pepmod_shift;
    vector[NpepmodCharges] pmc_shift;
    vector[Nmsruns*Npepmods] msrunXpepmod_shift;
    vector[Nmsruns*Npepmods] sampleXpepmod_shift;

    o_labu = rep_vector(base_labu, Nobservations);
    if (Nmsruns > 0) {
        msrun_shift_scale = msrun_shift_scale_tau * msrun_shift_scale_lambda_a / sqrt(msrun_shift_scale_lambda_t);
        msrun_shift_unscaled = msrun_eff2shift * msrun_effect;
        msrun_shift = msrun_shift_scale * msrun_shift_unscaled;
        o_labu = o_labu + msrun_shift[obs2msrun];
    }
    if (Nsamples > 0) {
        sample_shift_scale = sample_shift_scale_tau * sample_shift_scale_lambda_a / sqrt(sample_shift_scale_lambda_t);
        sample_shift_unscaled = sample_eff2shift * sample_effect;
        sample_shift = sample_shift_scale * sample_shift_unscaled;
        o_labu = o_labu + sample_shift[obs2sample];
    }
    if (Npepmods > 0) {
        pepmod_shift_scale = pepmod_shift_scale_tau * pepmod_shift_scale_lambda_a / sqrt(pepmod_shift_scale_lambda_t);
        pepmod_shift_unscaled = pepmod_eff2shift * pepmod_effect;
        pepmod_shift = pepmod_shift_scale * pepmod_shift_unscaled;
        o_labu = o_labu + pepmod_shift[obs2pepmod];
    }
    if (NpepmodCharges > 0) {
        pmc_shift_scale = pmc_shift_scale_tau * pmc_shift_scale_lambda_a / sqrt(pmc_shift_scale_lambda_t);
        pmc_shift_unscaled = pmc_eff2shift * pmc_effect;
        pmc_shift = pmc_shift_scale * pmc_shift_unscaled;
        o_labu = o_labu + pmc_shift[obs2pmc];
    }
    if (NmsrunXpepmodEffects > 0) {
        msrunXpepmod_shift_scale = msrunXpepmod_scale_tau * msrunXpepmod_shift_scale_lambda_a ./ sqrt(msrunXpepmod_shift_scale_lambda_t);
        msrunXpepmod_shift_unscaled = msrunXpepmod_eff2shift * msrunXpepmod_effect;
        msrunXpepmod_shift = msrunXpepmod_shift_scale .* msrunXpepmod_shift_unscaled;
        o_labu = o_labu + msrunXpepmod_shift[obs2msrunXpepmod];
    }
    if (NsampleXpepmodEffects > 0) {
        sampleXpepmod_shift_scale = sampleXpepmod_scale_tau * sampleXpepmod_shift_scale_lambda_a ./ sqrt(msrunXpepmod_shift_scale_lambda_t);
        sampleXpepmod_shift_unscaled = sampleXpepmod_eff2shift * sampleXpepmod_effect;
        sampleXpepmod_shift = sampleXpepmod_shift_scale .* sampleXpepmod_shift_unscaled;
        o_labu = o_labu + sampleXpepmod_shift[obs2sampleXpepmod];
    }
}

// msrun + msrun:experiment + pepmod + pepmod:charge + msrun:pepmod + experiment:pepmod
model {
    //base_labu ~ normal(0.0, 10.0);
    if (Npepmods > 0) {
        pepmod_shift_unscaled ~ normal(0.0, 1.0);
        pepmod_shift_scale_lambda_t ~ chi_square(pepmod_shift_scale_chi);
        pepmod_shift_scale_lambda_a ~ normal(0.0, 1.0); // 1.0 = 2/2
    }
    if (NpepmodCharges > 0) {
        pmc_shift_unscaled ~ normal(0.0, 1.0);
        pmc_shift_scale_lambda_t ~ chi_square(pmc_shift_scale_chi);
        pmc_shift_scale_lambda_a ~ normal(0.0, 1.0); // 1.0 = 2/2
    }
    if (Nmsruns > 0) {
        msrun_shift_unscaled ~ normal(0.0, 1.0);
        msrun_shift_scale_lambda_t ~ chi_square(msrun_shift_scale_chi);
        msrun_shift_scale_lambda_a ~ normal(0.0, 1.0); // 1.0 = 2/2
    }
    if (Nsamples > 0) {
        sample_shift_unscaled ~ normal(0.0, 1.0);
        sample_shift_scale_lambda_t ~ chi_square(sample_shift_scale_chi);
        sample_shift_scale_lambda_a ~ normal(0.0, 1.0); // 1.0 = 2/2
    }
    if (NmsrunXpepmodEffects > 0) {
      msrunXpepmod_shift_unscaled ~ normal(0.0, 1.0);
      msrunXpepmod_shift_scale_lambda_a ~ normal(0.0, 1.0); // 1.0 = 2/2
      msrunXpepmod_shift_scale_lambda_t ~ chi_square(msrunXpepmod_scale_chi);
    }
    if (NsampleXpepmodEffects > 0) {
      sampleXpepmod_shift_unscaled ~ normal(0.0, 1.0);
      sampleXpepmod_shift_scale_lambda_a ~ normal(0.0, 1.0); // 1.0 = 2/2
      sampleXpepmod_shift_scale_lambda_t ~ chi_square(sampleXpepmod_scale_chi);
    }

    // calculate the likelihood
    {
        vector[Nquanted] q_labu;
        vector[Nmissed] m_labu;

        q_labu = o_labu[quant2obs] + qShifts;
        m_labu = o_labu[miss2obs] + mShifts;

        qDataNorm ~ double_exponential(exp(q_labu - qLogStd), 1);
        // model quantitations and missing data
        1 ~ bernoulli_logit(q_labu * (zScale*zDetectionFactor) + (-mzShift * zScale*zDetectionFactor + zDetectionIntercept));
        0 ~ bernoulli_logit(m_labu * (zScale*zDetectionFactor) + (-mzShift * zScale*zDetectionFactor + zDetectionIntercept));
    }
}
