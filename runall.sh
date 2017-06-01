#!/usr/bin/env bash
# Script running all plotting macros.

# Closure tests
# ------------------------------------------------------------------------------

if false; then
    echo "RUNNING CLOSURE TESTS"

    python examples/closure.py --mass   0 --window 0.0 --save > logs/log_closure_0GeV_pm0.txt
    
    python examples/closure.py --mass  80 --window 0.2 --save > logs/log_closure_80GeV_pm20.txt
    python examples/closure.py --mass 100 --window 0.2 --save > logs/log_closure_100GeV_pm20.txt
    python examples/closure.py --mass 130 --window 0.2 --save > logs/log_closure_130GeV_pm20.txt
    python examples/closure.py --mass 160 --window 0.2 --save > logs/log_closure_160GeV_pm20.txt
    python examples/closure.py --mass 190 --window 0.2 --save > logs/log_closure_190GeV_pm20.txt
    python examples/closure.py --mass 220 --window 0.2 --save > logs/log_closure_220GeV_pm20.txt
    
    python examples/closure.py --mass  80 --window 0.3 --save > logs/log_closure_80GeV_pm30.txt
    python examples/closure.py --mass 100 --window 0.3 --save > logs/log_closure_100GeV_pm30.txt
    python examples/closure.py --mass 130 --window 0.3 --save > logs/log_closure_130GeV_pm30.txt
    python examples/closure.py --mass 160 --window 0.3 --save > logs/log_closure_160GeV_pm30.txt
    python examples/closure.py --mass 190 --window 0.3 --save > logs/log_closure_190GeV_pm30.txt
    python examples/closure.py --mass 220 --window 0.3 --save > logs/log_closure_220GeV_pm30.txt
fi


# Signal injection tests
# ------------------------------------------------------------------------------

if false; then
    echo "RUNNING SIGNAL INJECTION TESTS"

    python examples/signalinjection.py --mass 100 --save --inject > logs/log_signalinjection_100GeV_pm20_injected.txt
    python examples/signalinjection.py --mass 130 --save --inject > logs/log_signalinjection_130GeV_pm20_injected.txt
    python examples/signalinjection.py --mass 160 --save --inject > logs/log_signalinjection_160GeV_pm20_injected.txt
    python examples/signalinjection.py --mass 190 --save --inject > logs/log_signalinjection_190GeV_pm20_injected.txt
    python examples/signalinjection.py --mass 220 --save --inject > logs/log_signalinjection_220GeV_pm20_injected.txt
    
    python examples/signalinjection.py --mass 100 --save --inject --toys > logs/log_signalinjection_toys_100GeV_pm20_injected.txt
    python examples/signalinjection.py --mass 130 --save --inject --toys > logs/log_signalinjection_toys_130GeV_pm20_injected.txt
    python examples/signalinjection.py --mass 160 --save --inject --toys > logs/log_signalinjection_toys_160GeV_pm20_injected.txt
    python examples/signalinjection.py --mass 190 --save --inject --toys > logs/log_signalinjection_toys_190GeV_pm20_injected.txt
    python examples/signalinjection.py --mass 220 --save --inject --toys > logs/log_signalinjection_toys_220GeV_pm20_injected.txt
    
    python examples/signalinjection.py --mass 100 --save > logs/log_signalinjection_100GeV_pm20_notinjected.txt
    python examples/signalinjection.py --mass 130 --save > logs/log_signalinjection_130GeV_pm20_notinjected.txt
    python examples/signalinjection.py --mass 160 --save > logs/log_signalinjection_160GeV_pm20_notinjected.txt
    python examples/signalinjection.py --mass 190 --save > logs/log_signalinjection_190GeV_pm20_notinjected.txt
    python examples/signalinjection.py --mass 220 --save > logs/log_signalinjection_220GeV_pm20_notinjected.txt
    
    python examples/signalinjection.py --mass 100 --save --toys > logs/log_signalinjection_toys_100GeV_pm20_notinjected.txt
    python examples/signalinjection.py --mass 130 --save --toys > logs/log_signalinjection_toys_130GeV_pm20_notinjected.txt
    python examples/signalinjection.py --mass 160 --save --toys > logs/log_signalinjection_toys_160GeV_pm20_notinjected.txt
    python examples/signalinjection.py --mass 190 --save --toys > logs/log_signalinjection_toys_190GeV_pm20_notinjected.txt
    python examples/signalinjection.py --mass 220 --save --toys > logs/log_signalinjection_toys_220GeV_pm20_notinjected.txt
fi


# Generating TF outputs
# ------------------------------------------------------------------------------

if true; then
    echo "GENERATING TF OUTPUTS"
    #python examples/tf.py --mass  80 --save                                 > logs/log_tf_80GeV_pm20.txt
    #python examples/tf.py --mass  85 --save                                 > logs/log_tf_85GeV_pm20.txt
    #python examples/tf.py --mass  90 --save                                 > logs/log_tf_90GeV_pm20.txt
    #python examples/tf.py --mass 100 --save --subtractWZMC --subtractWZdata > logs/log_tf_100GeV_pm20_subWZ.txt
    #python examples/tf.py --mass 130 --save --subtractWZMC --subtractWZdata > logs/log_tf_130GeV_pm20_subWZ.txt
    #python examples/tf.py --mass 160 --save --subtractWZMC --subtractWZdata > logs/log_tf_160GeV_pm20_subWZ.txt
    #python examples/tf.py --mass 190 --save --subtractWZMC --subtractWZdata > logs/log_tf_190GeV_pm20_subWZ.txt
    #python examples/tf.py --mass 220 --save --subtractWZMC --subtractWZdata > logs/log_tf_220GeV_pm20_subWZ.txt
    
    python examples/tf.py --mass 100 --save > logs/log_tf_100GeV_pm20.txt
    python examples/tf.py --mass 130 --save > logs/log_tf_130GeV_pm20.txt
    python examples/tf.py --mass 160 --save > logs/log_tf_160GeV_pm20.txt
    python examples/tf.py --mass 190 --save > logs/log_tf_190GeV_pm20.txt
    python examples/tf.py --mass 220 --save > logs/log_tf_220GeV_pm20.txt
    
    #python examples/tf.py --mass  80 --save --window 0.3                                 > logs/log_tf_80GeV_pm30.txt
    #python examples/tf.py --mass  85 --save --window 0.3                                 > logs/log_tf_85GeV_pm30.txt
    #python examples/tf.py --mass  90 --save --window 0.3                                 > logs/log_tf_90GeV_pm30.txt
    #python examples/tf.py --mass 100 --save --window 0.3 --subtractWZMC --subtractWZdata > logs/log_tf_100GeV_pm30_subWZ.txt
    #python examples/tf.py --mass 130 --save --window 0.3 --subtractWZMC --subtractWZdata > logs/log_tf_130GeV_pm30_subWZ.txt
    #python examples/tf.py --mass 160 --save --window 0.3 --subtractWZMC --subtractWZdata > logs/log_tf_160GeV_pm30_subWZ.txt
    #python examples/tf.py --mass 190 --save --window 0.3 --subtractWZMC --subtractWZdata > logs/log_tf_190GeV_pm30_subWZ.txt
    #python examples/tf.py --mass 220 --save --window 0.3 --subtractWZMC --subtractWZdata > logs/log_tf_220GeV_pm30_subWZ.txt
fi


# W/Z search
# ------------------------------------------------------------------------------

if false; then
    echo "RUNNING W/Z SEARCH"

    python examples/wzsearch.py --save              > logs/log_wzsearch_pm20.txt

    python examples/wzsearch.py --save --window 0.3 > logs/log_wzsearch_pm30.txt
fi


# Global background shape
# ------------------------------------------------------------------------------

if false; then
    echo "RUNNING GLOBAL BACKGROUND SHAPE"

    python examples/globalbackground.py --mass 100 --save --inject > logs/log_globalbackground_100GeV_injected.txt
    python examples/globalbackground.py --mass 130 --save --inject > logs/log_globalbackground_130GeV_injected.txt
    python examples/globalbackground.py --mass 160 --save --inject > logs/log_globalbackground_160GeV_injected.txt
    python examples/globalbackground.py --mass 190 --save --inject > logs/log_globalbackground_190GeV_injected.txt
    python examples/globalbackground.py --mass 220 --save --inject > logs/log_globalbackground_220GeV_injected.txt
    
    python examples/globalbackground.py --mass 100 --save          > logs/log_globalbackground_100GeV_notinjected.txt
    python examples/globalbackground.py --mass 130 --save          > logs/log_globalbackground_130GeV_notinjected.txt
    python examples/globalbackground.py --mass 160 --save          > logs/log_globalbackground_160GeV_notinjected.txt
    python examples/globalbackground.py --mass 190 --save          > logs/log_globalbackground_190GeV_notinjected.txt
    python examples/globalbackground.py --mass 220 --save          > logs/log_globalbackground_220GeV_notinjected.txt
fi
