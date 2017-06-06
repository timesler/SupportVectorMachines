/*
################################################################################
## HPCC SYSTEMS software Copyright (C) 2016 HPCC Systems®.  All rights reserved.
################################################################################
*/
/**
  * Data and model preparation for performance testing of SupportVectorMachines
  * bundle.
  */

IMPORT PBblas;
IMPORT PBblas.internal as int;
IMPORT int.Types as iTypes;
IMPORT PBblas.Types;
IMPORT int.MatDims;
IMPORT PBblas.test as Tests;
IMPORT Tests.MakeTestMatrix as tm;
IMPORT ML_Core as ML;
IMPORT ML.Types as Core_Types;
IMPORT SupportVectorMachines as SVM;

Layout_Cell := Types.Layout_Cell;
NumericField := Core_Types.NumericField;
Layout_Model := Core_Types.Layout_Model;

// Test configuration Parameters -- modify these to vary the testing

N := 1000;   // Number of rows in X
M := 20;     // Number of columns in X

// Small medium and large sizes
Myriad_Count_1    := 1;
Myriad_Count_2    := 2;
Myriad_Count_4    := 4;
Myriad_Count_8    := 8;
Myriad_Count_16   := 16;
Myriad_Count_20   := 20;
Myriad_Count_21   := 21;
Myriad_Count_22   := 22;
Myriad_Count_64   := 64;
Myriad_Count_256  := 256;
Myriad_Count_1024 := 1024;

density := 1.0; // 1.0 is fully dense. 0.0 is empty.
// End of config parameters

// Generate test data for X and Y

X_LC := tm.RandomMatrix(N, M, density, 1);
Y_LC := tm.RandomMatrix(N, 1, density, 1);
Y_min := MIN(Y_LC, v);
Y_max := MAX(Y_LC, v);

NumericField asNumericFieldX(Layout_Cell L) := TRANSFORM
  SELF.wi     := L.wi_id;
  SELF.id     := L.x;
  SELF.number := L.y;
  SELF.value  := L.v;
END;
NumericField asNumericFieldY(Layout_Cell L) := TRANSFORM
  SELF.wi     := L.wi_id;
  SELF.id     := L.x;
  SELF.number := L.y;
  SELF.value  := ROUND((L.v - Y_min) / (Y_max - Y_min)) * 2 - 1;
END;

X := PROJECT(X_LC, asNumericFieldX(LEFT));
Y := PROJECT(Y_LC, asNumericFieldY(LEFT));

NumericField make_myriad(NumericField d, UNSIGNED c) := TRANSFORM
  SELF.wi := c;
  SELF    := d;
END;

X_myr_1    := NORMALIZE(X, Myriad_Count_1, make_myriad(LEFT, COUNTER));
X_myr_2    := NORMALIZE(X, Myriad_Count_2, make_myriad(LEFT, COUNTER));
X_myr_4    := NORMALIZE(X, Myriad_Count_4, make_myriad(LEFT, COUNTER));
X_myr_8    := NORMALIZE(X, Myriad_Count_8, make_myriad(LEFT, COUNTER));
X_myr_16   := NORMALIZE(X, Myriad_Count_16, make_myriad(LEFT, COUNTER));
X_myr_20   := NORMALIZE(X, Myriad_Count_20, make_myriad(LEFT, COUNTER));
X_myr_21   := NORMALIZE(X, Myriad_Count_21, make_myriad(LEFT, COUNTER));
X_myr_22   := NORMALIZE(X, Myriad_Count_22, make_myriad(LEFT, COUNTER));
X_myr_64   := NORMALIZE(X, Myriad_Count_64, make_myriad(LEFT, COUNTER));
X_myr_256  := NORMALIZE(X, Myriad_Count_256, make_myriad(LEFT, COUNTER));
X_myr_1024 := NORMALIZE(X, Myriad_Count_1024, make_myriad(LEFT, COUNTER));
Y_myr_1    := NORMALIZE(Y, Myriad_Count_1, make_myriad(LEFT, COUNTER));
Y_myr_2    := NORMALIZE(Y, Myriad_Count_2, make_myriad(LEFT, COUNTER));
Y_myr_4    := NORMALIZE(Y, Myriad_Count_4, make_myriad(LEFT, COUNTER));
Y_myr_8    := NORMALIZE(Y, Myriad_Count_8, make_myriad(LEFT, COUNTER));
Y_myr_16   := NORMALIZE(Y, Myriad_Count_16, make_myriad(LEFT, COUNTER));
Y_myr_20   := NORMALIZE(Y, Myriad_Count_20, make_myriad(LEFT, COUNTER));
Y_myr_21   := NORMALIZE(Y, Myriad_Count_21, make_myriad(LEFT, COUNTER));
Y_myr_22   := NORMALIZE(Y, Myriad_Count_22, make_myriad(LEFT, COUNTER));
Y_myr_64   := NORMALIZE(Y, Myriad_Count_64, make_myriad(LEFT, COUNTER));
Y_myr_256  := NORMALIZE(Y, Myriad_Count_256, make_myriad(LEFT, COUNTER));
Y_myr_1024 := NORMALIZE(Y, Myriad_Count_1024, make_myriad(LEFT, COUNTER));

OUTPUT(X_myr_1,,'GetModelPerfMyr_X_1.dat', OVERWRITE);
OUTPUT(X_myr_2,,'GetModelPerfMyr_X_2.dat', OVERWRITE);
OUTPUT(X_myr_4,,'GetModelPerfMyr_X_4.dat', OVERWRITE);
OUTPUT(X_myr_8,,'GetModelPerfMyr_X_8.dat', OVERWRITE);
OUTPUT(X_myr_16,,'GetModelPerfMyr_X_16.dat', OVERWRITE);
OUTPUT(X_myr_20,,'GetModelPerfMyr_X_20.dat', OVERWRITE);
OUTPUT(X_myr_21,,'GetModelPerfMyr_X_21.dat', OVERWRITE);
OUTPUT(X_myr_22,,'GetModelPerfMyr_X_22.dat', OVERWRITE);
OUTPUT(X_myr_64,,'GetModelPerfMyr_X_64.dat', OVERWRITE);
OUTPUT(X_myr_256,,'GetModelPerfMyr_X_256.dat', OVERWRITE);
OUTPUT(X_myr_1024,,'GetModelPerfMyr_X_1024.dat', OVERWRITE);
OUTPUT(Y_myr_1,,'GetModelPerfMyr_Y_1.dat', OVERWRITE);
OUTPUT(Y_myr_2,,'GetModelPerfMyr_Y_2.dat', OVERWRITE);
OUTPUT(Y_myr_4,,'GetModelPerfMyr_Y_4.dat', OVERWRITE);
OUTPUT(Y_myr_8,,'GetModelPerfMyr_Y_8.dat', OVERWRITE);
OUTPUT(Y_myr_16,,'GetModelPerfMyr_Y_16.dat', OVERWRITE);
OUTPUT(Y_myr_20,,'GetModelPerfMyr_Y_20.dat', OVERWRITE);
OUTPUT(Y_myr_21,,'GetModelPerfMyr_Y_21.dat', OVERWRITE);
OUTPUT(Y_myr_22,,'GetModelPerfMyr_Y_22.dat', OVERWRITE);
OUTPUT(Y_myr_64,,'GetModelPerfMyr_Y_64.dat', OVERWRITE);
OUTPUT(Y_myr_256,,'GetModelPerfMyr_Y_256.dat', OVERWRITE);
OUTPUT(Y_myr_1024,,'GetModelPerfMyr_Y_1024.dat', OVERWRITE);

SVCModel_1    := SVM.SVC().GetModel(X_myr_1, ML.Discretize.ByRounding(Y_myr_1));
SVCModel_2    := SVM.SVC().GetModel(X_myr_2, ML.Discretize.ByRounding(Y_myr_2));
SVCModel_4    := SVM.SVC().GetModel(X_myr_4, ML.Discretize.ByRounding(Y_myr_4));
SVCModel_8    := SVM.SVC().GetModel(X_myr_8, ML.Discretize.ByRounding(Y_myr_8));
SVCModel_16   := SVM.SVC().GetModel(X_myr_16, ML.Discretize.ByRounding(Y_myr_16));
SVCModel_20   := SVM.SVC().GetModel(X_myr_20, ML.Discretize.ByRounding(Y_myr_20));
SVCModel_21   := SVM.SVC().GetModel(X_myr_21, ML.Discretize.ByRounding(Y_myr_21));
SVCModel_22   := SVM.SVC().GetModel(X_myr_22, ML.Discretize.ByRounding(Y_myr_22));
SVCModel_64   := SVM.SVC().GetModel(X_myr_64, ML.Discretize.ByRounding(Y_myr_64));
SVCModel_256  := SVM.SVC().GetModel(X_myr_256, ML.Discretize.ByRounding(Y_myr_256));
SVCModel_1024 := SVM.SVC().GetModel(X_myr_1024, ML.Discretize.ByRounding(Y_myr_1024));

SVRModel_1 := SVM.SVR(NAMED X := X_myr_1, NAMED Y := Y_myr_1).GetModel();
SVRModel_2 := SVM.SVR(NAMED X := X_myr_2, NAMED Y := Y_myr_2).GetModel();
SVRModel_4 := SVM.SVR(NAMED X := X_myr_4, NAMED Y := Y_myr_4).GetModel();
SVRModel_8 := SVM.SVR(NAMED X := X_myr_8, NAMED Y := Y_myr_8).GetModel();
SVRModel_16 := SVM.SVR(NAMED X := X_myr_16, NAMED Y := Y_myr_16).GetModel();
SVRModel_20 := SVM.SVR(NAMED X := X_myr_20, NAMED Y := Y_myr_20).GetModel();
SVRModel_21 := SVM.SVR(NAMED X := X_myr_21, NAMED Y := Y_myr_21).GetModel();
SVRModel_22 := SVM.SVR(NAMED X := X_myr_22, NAMED Y := Y_myr_22).GetModel();
SVRModel_64 := SVM.SVR(NAMED X := X_myr_64, NAMED Y := Y_myr_64).GetModel();
SVRModel_256 := SVM.SVR(NAMED X := X_myr_256, NAMED Y := Y_myr_256).GetModel();
SVRModel_1024 := SVM.SVR(NAMED X := X_myr_1024, NAMED Y := Y_myr_1024).GetModel();

OUTPUT(SVCModel_1,,'PredictPerfMyr_SVC_1.dat', OVERWRITE);
OUTPUT(SVCModel_2,,'PredictPerfMyr_SVC_2.dat', OVERWRITE);
OUTPUT(SVCModel_4,,'PredictPerfMyr_SVC_4.dat', OVERWRITE);
OUTPUT(SVCModel_8,,'PredictPerfMyr_SVC_8.dat', OVERWRITE);
OUTPUT(SVCModel_16,,'PredictPerfMyr_SVC_16.dat', OVERWRITE);
OUTPUT(SVCModel_20,,'PredictPerfMyr_SVC_20.dat', OVERWRITE);
OUTPUT(SVCModel_21,,'PredictPerfMyr_SVC_21.dat', OVERWRITE);
OUTPUT(SVCModel_22,,'PredictPerfMyr_SVC_22.dat', OVERWRITE);
OUTPUT(SVCModel_64,,'PredictPerfMyr_SVC_64.dat', OVERWRITE);
OUTPUT(SVCModel_256,,'PredictPerfMyr_SVC_256.dat', OVERWRITE);
OUTPUT(SVCModel_1024,,'PredictPerfMyr_SVC_1024.dat', OVERWRITE);
OUTPUT(SVRModel_1,,'PredictPerfMyr_SVR_1.dat', OVERWRITE);
OUTPUT(SVRModel_2,,'PredictPerfMyr_SVR_2.dat', OVERWRITE);
OUTPUT(SVRModel_4,,'PredictPerfMyr_SVR_4.dat', OVERWRITE);
OUTPUT(SVRModel_8,,'PredictPerfMyr_SVR_8.dat', OVERWRITE);
OUTPUT(SVRModel_16,,'PredictPerfMyr_SVR_16.dat', OVERWRITE);
OUTPUT(SVRModel_20,,'PredictPerfMyr_SVR_20.dat', OVERWRITE);
OUTPUT(SVRModel_21,,'PredictPerfMyr_SVR_21.dat', OVERWRITE);
OUTPUT(SVRModel_22,,'PredictPerfMyr_SVR_22.dat', OVERWRITE);
OUTPUT(SVRModel_64,,'PredictPerfMyr_SVR_64.dat', OVERWRITE);
OUTPUT(SVRModel_256,,'PredictPerfMyr_SVR_256.dat', OVERWRITE);
OUTPUT(SVRModel_1024,,'PredictPerfMyr_SVR_1024.dat', OVERWRITE);
