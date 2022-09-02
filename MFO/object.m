function f=object(x)
load modeldata01
load modeldata02
load modeldata_c
x=[model_Q x];
[Predict_1 MSE1] = predictor(x, dmodel01);
[Predict_2 MSE2] = predictor(x, dmodel02);
f(1)=abs(Predict_1-model_P);
f(2)=-Predict_2;
end
