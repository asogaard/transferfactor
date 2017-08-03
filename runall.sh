#!/usr/bin/env bash
# Script running all plotting macros.

# Flags
DECORRELATION=false
CLOSURE=false
INJECTION=false
TF_OUTPUTS=false
WZ_SEARCH=false
GBS=false
GLOBALBACKGROUNDPLOTS=false
DATA_DISTRIBUTIONS=false
INTERPOLATION=true
ACCEPTANCE=false
CUTOPTIMISATION=false
PAPERFIGS=false
VALIDATION=false


# Substructure de-correlation
# ------------------------------------------------------------------------------

if $DECORRELATION; then
    echo "RUNNING SUBSTRUCTURE DE-CORRELATION"

    python examples/decorrelation.py --save  > logs/log_decorrelation.txt &
    wait
fi


# Closure tests
# ------------------------------------------------------------------------------

if $CLOSURE; then
    echo "RUNNING CLOSURE TESTS"

    python examples/closure.py --mass   0 --window 0.0 --save > logs/log_closure_0GeV_pm0.txt    &
    
    python examples/closure.py --mass  80 --window 0.2 --save > logs/log_closure_80GeV_pm20.txt  &
    python examples/closure.py --mass 100 --window 0.2 --save > logs/log_closure_100GeV_pm20.txt &
    python examples/closure.py --mass 130 --window 0.2 --save > logs/log_closure_130GeV_pm20.txt &
    python examples/closure.py --mass 160 --window 0.2 --save > logs/log_closure_160GeV_pm20.txt &
    python examples/closure.py --mass 190 --window 0.2 --save > logs/log_closure_190GeV_pm20.txt &
    python examples/closure.py --mass 220 --window 0.2 --save > logs/log_closure_220GeV_pm20.txt &
    wait

    python examples/closure.py --mass  80 --window 0.3 --save > logs/log_closure_80GeV_pm30.txt  &
    python examples/closure.py --mass 100 --window 0.3 --save > logs/log_closure_100GeV_pm30.txt &
    python examples/closure.py --mass 130 --window 0.3 --save > logs/log_closure_130GeV_pm30.txt &
    python examples/closure.py --mass 160 --window 0.3 --save > logs/log_closure_160GeV_pm30.txt &
    python examples/closure.py --mass 190 --window 0.3 --save > logs/log_closure_190GeV_pm30.txt &
    python examples/closure.py --mass 220 --window 0.3 --save > logs/log_closure_220GeV_pm30.txt &
    wait
fi


# Signal injection tests
# ------------------------------------------------------------------------------

if $INJECTION; then
    echo "RUNNING SIGNAL INJECTION TESTS"

    python examples/signalinjection.py --mass 100 --save --inject > logs/log_signalinjection_100GeV_pm20_injected.txt &
    python examples/signalinjection.py --mass 130 --save --inject > logs/log_signalinjection_130GeV_pm20_injected.txt &
    python examples/signalinjection.py --mass 160 --save --inject > logs/log_signalinjection_160GeV_pm20_injected.txt &
    python examples/signalinjection.py --mass 190 --save --inject > logs/log_signalinjection_190GeV_pm20_injected.txt &
    python examples/signalinjection.py --mass 220 --save --inject > logs/log_signalinjection_220GeV_pm20_injected.txt &
    wait
    
    python examples/signalinjection.py --mass 100 --save --inject --toys > logs/log_signalinjection_toys_100GeV_pm20_injected.txt &
    python examples/signalinjection.py --mass 130 --save --inject --toys > logs/log_signalinjection_toys_130GeV_pm20_injected.txt &
    python examples/signalinjection.py --mass 160 --save --inject --toys > logs/log_signalinjection_toys_160GeV_pm20_injected.txt &
    python examples/signalinjection.py --mass 190 --save --inject --toys > logs/log_signalinjection_toys_190GeV_pm20_injected.txt &
    python examples/signalinjection.py --mass 220 --save --inject --toys > logs/log_signalinjection_toys_220GeV_pm20_injected.txt &
    wait

    python examples/signalinjection.py --mass 100 --save > logs/log_signalinjection_100GeV_pm20_notinjected.txt &
    python examples/signalinjection.py --mass 130 --save > logs/log_signalinjection_130GeV_pm20_notinjected.txt &
    python examples/signalinjection.py --mass 160 --save > logs/log_signalinjection_160GeV_pm20_notinjected.txt &
    python examples/signalinjection.py --mass 190 --save > logs/log_signalinjection_190GeV_pm20_notinjected.txt &
    python examples/signalinjection.py --mass 220 --save > logs/log_signalinjection_220GeV_pm20_notinjected.txt &
    wait

    python examples/signalinjection.py --mass 100 --save --toys > logs/log_signalinjection_toys_100GeV_pm20_notinjected.txt &
    python examples/signalinjection.py --mass 130 --save --toys > logs/log_signalinjection_toys_130GeV_pm20_notinjected.txt &
    python examples/signalinjection.py --mass 160 --save --toys > logs/log_signalinjection_toys_160GeV_pm20_notinjected.txt &
    python examples/signalinjection.py --mass 190 --save --toys > logs/log_signalinjection_toys_190GeV_pm20_notinjected.txt &
    python examples/signalinjection.py --mass 220 --save --toys > logs/log_signalinjection_toys_220GeV_pm20_notinjected.txt &
    wait
fi


# Generating TF outputs
# ------------------------------------------------------------------------------

if $TF_OUTPUTS; then
    echo "GENERATING TF OUTPUTS"
    python examples/tf.py --mass   0 --save  --subtractWZMC --subtractWZdata > logs/log_tf_0GeV_pm20.txt   & sleep 5

    python examples/tf.py --mass  85              --save                     > logs/log_tf_85GeV_pm20.txt  & sleep 5
    python examples/tf.py --mass  85 --window 0.3 --save                     > logs/log_tf_85GeV_pm30.txt  & sleep 5

    python examples/tf.py --mass 100 --save  --subtractWZMC --subtractWZdata > logs/log_tf_100GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 105 --save  --subtractWZMC --subtractWZdata > logs/log_tf_105GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 110 --save  --subtractWZMC --subtractWZdata > logs/log_tf_110GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 115 --save  --subtractWZMC --subtractWZdata > logs/log_tf_115GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 120 --save  --subtractWZMC --subtractWZdata > logs/log_tf_120GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 125 --save  --subtractWZMC --subtractWZdata > logs/log_tf_125GeV_pm20.txt & sleep 5
    wait

    python examples/tf.py --mass 130 --save  --subtractWZMC --subtractWZdata > logs/log_tf_130GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 135 --save  --subtractWZMC --subtractWZdata > logs/log_tf_135GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 140 --save  --subtractWZMC --subtractWZdata > logs/log_tf_140GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 145 --save  --subtractWZMC --subtractWZdata > logs/log_tf_145GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 150 --save  --subtractWZMC --subtractWZdata > logs/log_tf_150GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 155 --save  --subtractWZMC --subtractWZdata > logs/log_tf_155GeV_pm20.txt & sleep 5
    wait

    python examples/tf.py --mass 160 --save  --subtractWZMC --subtractWZdata > logs/log_tf_160GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 165 --save  --subtractWZMC --subtractWZdata > logs/log_tf_165GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 170 --save  --subtractWZMC --subtractWZdata > logs/log_tf_170GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 175 --save  --subtractWZMC --subtractWZdata > logs/log_tf_175GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 180 --save  --subtractWZMC --subtractWZdata > logs/log_tf_180GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 185 --save  --subtractWZMC --subtractWZdata > logs/log_tf_185GeV_pm20.txt & sleep 5
    wait

    python examples/tf.py --mass 190 --save  --subtractWZMC --subtractWZdata > logs/log_tf_190GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 195 --save  --subtractWZMC --subtractWZdata > logs/log_tf_195GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 200 --save  --subtractWZMC --subtractWZdata > logs/log_tf_200GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 205 --save  --subtractWZMC --subtractWZdata > logs/log_tf_205GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 210 --save  --subtractWZMC --subtractWZdata > logs/log_tf_210GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 215 --save  --subtractWZMC --subtractWZdata > logs/log_tf_215GeV_pm20.txt & sleep 5
    wait

    python examples/tf.py --mass 220 --save  --subtractWZMC --subtractWZdata > logs/log_tf_220GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 225 --save  --subtractWZMC --subtractWZdata > logs/log_tf_225GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 230 --save  --subtractWZMC --subtractWZdata > logs/log_tf_230GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 235 --save  --subtractWZMC --subtractWZdata > logs/log_tf_235GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 240 --save  --subtractWZMC --subtractWZdata > logs/log_tf_240GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 245 --save  --subtractWZMC --subtractWZdata > logs/log_tf_245GeV_pm20.txt & sleep 5
    python examples/tf.py --mass 250 --save  --subtractWZMC --subtractWZdata > logs/log_tf_250GeV_pm20.txt & sleep 5
    wait
fi


# W/Z search
# ------------------------------------------------------------------------------

if $WZ_SEARCH; then
    echo "RUNNING W/Z SEARCH"

    python examples/wzsearch.py --save              > logs/log_wzsearch_pm20.txt &

    python examples/wzsearch.py --save --window 0.3 > logs/log_wzsearch_pm30.txt &
    wait
fi


# Global background shape - production
# ------------------------------------------------------------------------------

if $GBS; then
    echo "RUNNING GLOBAL BACKGROUND SHAPE PRODUCTION"

    python examples/gbs.py --data --subtractWZMC --subtractWZdata --save  > logs/log_gbs_data_subtractWZMC_subtractMCdata.txt &
    python examples/gbs.py                                        --save  > logs/log_gbs_MC_subtractWZMC_subtractMCdata.txt   &
    wait
fi


# Global background shape - plotting
# ------------------------------------------------------------------------------

if $GLOBALBACKGROUNDPLOTS; then
    echo "RUNNING GLOBAL BACKGROUND SHAPE PLOTTING"

    #python examples/globalbackground.py --mass 100 --save --data          > logs/log_globalbackground_100GeV_data_notinjected.txt &
    #python examples/globalbackground.py --mass 130 --save --data          > logs/log_globalbackground_130GeV_data_notinjected.txt &
    #python examples/globalbackground.py --mass 160 --save --data          > logs/log_globalbackground_160GeV_data_notinjected.txt &
    #python examples/globalbackground.py --mass 190 --save --data          > logs/log_globalbackground_190GeV_data_notinjected.txt &
    #python examples/globalbackground.py --mass 220 --save --data          > logs/log_globalbackground_220GeV_data_notinjected.txt &

    #python examples/globalbackground.py --mass 100 --save --data --inject > logs/log_globalbackground_100GeV_data_injected.txt &
    #python examples/globalbackground.py --mass 130 --save --data --inject > logs/log_globalbackground_130GeV_data_injected.txt &
    #python examples/globalbackground.py --mass 160 --save --data --inject > logs/log_globalbackground_160GeV_data_injected.txt &
    #python examples/globalbackground.py --mass 190 --save --data --inject > logs/log_globalbackground_190GeV_data_injected.txt &
    #python examples/globalbackground.py --mass 220 --save --data --inject > logs/log_globalbackground_220GeV_data_injected.txt &
    wait

    python examples/globalbackground.py --mass 100 --save --inject > logs/log_globalbackground_100GeV_injected.txt &
    python examples/globalbackground.py --mass 130 --save --inject > logs/log_globalbackground_130GeV_injected.txt &
    python examples/globalbackground.py --mass 160 --save --inject > logs/log_globalbackground_160GeV_injected.txt &
    python examples/globalbackground.py --mass 190 --save --inject > logs/log_globalbackground_190GeV_injected.txt &
    python examples/globalbackground.py --mass 220 --save --inject > logs/log_globalbackground_220GeV_injected.txt &
    
    python examples/globalbackground.py --mass 100 --save          > logs/log_globalbackground_100GeV_notinjected.txt &
    python examples/globalbackground.py --mass 130 --save          > logs/log_globalbackground_130GeV_notinjected.txt &
    python examples/globalbackground.py --mass 160 --save          > logs/log_globalbackground_160GeV_notinjected.txt &
    python examples/globalbackground.py --mass 190 --save          > logs/log_globalbackground_190GeV_notinjected.txt &
    python examples/globalbackground.py --mass 220 --save          > logs/log_globalbackground_220GeV_notinjected.txt &
    wait
fi


# Unblinded data distributions
# ------------------------------------------------------------------------------

if $DATA_DISTRIBUTIONS; then
    echo "RUNNING UNBLINDED DATA DISTRIBUTIONS"

    python examples/datadistributions.py --mass   0 --save > logs/log_datadistributions_0GeV.txt &

    python examples/datadistributions.py --mass 100 --save > logs/log_datadistributions_100GeV.txt &
    python examples/datadistributions.py --mass 105 --save > logs/log_datadistributions_105GeV.txt &
    python examples/datadistributions.py --mass 110 --save > logs/log_datadistributions_110GeV.txt &
    python examples/datadistributions.py --mass 115 --save > logs/log_datadistributions_115GeV.txt &
    python examples/datadistributions.py --mass 120 --save > logs/log_datadistributions_120GeV.txt &
    python examples/datadistributions.py --mass 125 --save > logs/log_datadistributions_125GeV.txt &
    wait

    python examples/datadistributions.py --mass 130 --save > logs/log_datadistributions_130GeV.txt &
    python examples/datadistributions.py --mass 135 --save > logs/log_datadistributions_135GeV.txt &
    python examples/datadistributions.py --mass 140 --save > logs/log_datadistributions_140GeV.txt &
    python examples/datadistributions.py --mass 145 --save > logs/log_datadistributions_145GeV.txt &
    python examples/datadistributions.py --mass 150 --save > logs/log_datadistributions_150GeV.txt &
    python examples/datadistributions.py --mass 155 --save > logs/log_datadistributions_155GeV.txt &
    wait

    python examples/datadistributions.py --mass 160 --save > logs/log_datadistributions_160GeV.txt &
    python examples/datadistributions.py --mass 165 --save > logs/log_datadistributions_165GeV.txt &
    python examples/datadistributions.py --mass 170 --save > logs/log_datadistributions_170GeV.txt &
    python examples/datadistributions.py --mass 175 --save > logs/log_datadistributions_175GeV.txt &
    python examples/datadistributions.py --mass 180 --save > logs/log_datadistributions_180GeV.txt &
    python examples/datadistributions.py --mass 185 --save > logs/log_datadistributions_185GeV.txt &
    wait

    python examples/datadistributions.py --mass 190 --save > logs/log_datadistributions_190GeV.txt &
    python examples/datadistributions.py --mass 195 --save > logs/log_datadistributions_195GeV.txt &
    python examples/datadistributions.py --mass 200 --save > logs/log_datadistributions_200GeV.txt &
    python examples/datadistributions.py --mass 205 --save > logs/log_datadistributions_205GeV.txt &
    python examples/datadistributions.py --mass 210 --save > logs/log_datadistributions_210GeV.txt &
    python examples/datadistributions.py --mass 215 --save > logs/log_datadistributions_215GeV.txt &
    wait

    python examples/datadistributions.py --mass 220 --save > logs/log_datadistributions_220GeV.txt &
    python examples/datadistributions.py --mass 225 --save > logs/log_datadistributions_225GeV.txt &
    python examples/datadistributions.py --mass 230 --save > logs/log_datadistributions_230GeV.txt &
    python examples/datadistributions.py --mass 235 --save > logs/log_datadistributions_235GeV.txt &
    python examples/datadistributions.py --mass 240 --save > logs/log_datadistributions_240GeV.txt &
    python examples/datadistributions.py --mass 245 --save > logs/log_datadistributions_245GeV.txt &
    python examples/datadistributions.py --mass 250 --save > logs/log_datadistributions_250GeV.txt &
    wait
fi


# Signal shape interpolation
# ------------------------------------------------------------------------------

if $INTERPOLATION; then
    echo "RUNNING SIGNAL INTERPOLATION"

    python examples/interpolation.py          --save > logs/log_interpolation_isrgamma.txt &
    python examples/interpolation.py --isrjet --save > logs/log_interpolation_isrjet.txt   &
    #wait
fi


# Signal acceptance
# ------------------------------------------------------------------------------

if $ACCEPTANCE; then
    echo "RUNNING SIGNAL ACCEPTANCE"

    python examples/acceptance.py --save > logs/log_acceptance.txt &
    wait
fi


# tau21DDT cut optimisation
# ------------------------------------------------------------------------------

if $CUTOPTIMISATION; then
    echo "RUNNING TAU21DDT CUT OPTIMISATION"

    python examples/cutoptimisation.py --save > logs/log_cutoptimisation.txt &
    wait
fi


# Paper plots
# ------------------------------------------------------------------------------

if $PAPERFIGS; then
    echo "RUNNING PAPER FIGURES"

    python examples/paperfigures.py --save > logs/log_paperfigures.txt &
    wait
fi


# Transfer factor fit validation
# ------------------------------------------------------------------------------

if $VALIDATION; then
    echo "RUNNING TF VALIDATION"

    python examples/validation.py --N 10 --save > logs/log_validation.txt &


    #python examples/validation.py --mass  85 --N 10 --save > logs/log_validation_85GeV.txt &
    #wait 

    #python examples/validation.py --mass 100 --N 10 --save > logs/log_validation_100GeV.txt &
    #python examples/validation.py --mass 110 --N 10 --save > logs/log_validation_110GeV.txt &
    #python examples/validation.py --mass 120 --N 10 --save > logs/log_validation_120GeV.txt &
    #python examples/validation.py --mass 130 --N 10 --save > logs/log_validation_130GeV.txt &
    #python examples/validation.py --mass 140 --N 10 --save > logs/log_validation_140GeV.txt &
    #python examples/validation.py --mass 150 --N 10 --save > logs/log_validation_150GeV.txt &
    #python examples/validation.py --mass 160 --N 10 --save > logs/log_validation_160GeV.txt &
    #python examples/validation.py --mass 170 --N 10 --save > logs/log_validation_170GeV.txt &
    #wait

    #python examples/validation.py --mass 180 --N 10 --save > logs/log_validation_180GeV.txt &
    #python examples/validation.py --mass 190 --N 10 --save > logs/log_validation_190GeV.txt &
    #python examples/validation.py --mass 200 --N 10 --save > logs/log_validation_200GeV.txt &
    #python examples/validation.py --mass 210 --N 10 --save > logs/log_validation_210GeV.txt &
    #python examples/validation.py --mass 220 --N 10 --save > logs/log_validation_220GeV.txt &
    #python examples/validation.py --mass 230 --N 10 --save > logs/log_validation_230GeV.txt &
    #python examples/validation.py --mass 240 --N 10 --save > logs/log_validation_240GeV.txt &
    #python examples/validation.py --mass 250 --N 10 --save > logs/log_validation_250GeV.txt &
    #wait
fi
