// Converted to and from SVM and generic formats for ML/Classify
IMPORT ML_Core as ML;
IMPORT $ as SVM;
IMPORT ML_Core.Types AS ML_Types;

/**
 * Module for various data conversions, including from a NumericField to an
 * SVM_Instance, and to/from a Layout_Model/NumericField.
 */
EXPORT Converted := MODULE
  //Instance data
  // Don't need from instance.  This can be a one-way trip.
  dummy := DATASET([], ML_Types.NumericField);

  /**
   * Convert dataset in NumericField format (with separate NumericFields for
   * dependent and independent variables) to an SVM_Instance format.
   * @param Ind NumericField dataset of independent variables.
   * @param Dep NumericField dataset of dependent variable(s) (default: empty dataset).
   * @return Dataset converted to SVM_Instance format (see SupportVectorMachines.Types
   * for type definition).
   */
  EXPORT ToInstance(DATASET(ML_Types.NumericField) Ind,
                    DATASET(ML_Types.NumericField) Dep=dummy) := FUNCTION
    g_i := GROUP(SORT(Ind, wi, id, number), wi, id);
    SVM.Types.SVM_Feature cvtF(ML_Types.NumericField nf) := TRANSFORM
      SELF.nominal := nf.number;
      SELF.v := nf.value;
    END;
    SVM.Types.SVM_Instance rollNF(ML_Types.NumericField f,
                                  DATASET(ML_Types.NumericField) ds) := TRANSFORM
      SELF.wi := f.wi;
      SELF.rid := f.id;
      SELF.y := 0;
      SELF.max_value := MAX(ds,value);
      SELF.x := PROJECT(ds, cvtF(LEFT));
    END;
    indy := ROLLUP(g_i, GROUP, rollNF(LEFT, ROWS(LEFT)));
    SVM.Types.SVM_Instance getD(SVM.Types.SVM_Instance i, ML_Types.NumericField d):=TRANSFORM
      SELF.y := d.value;
      SELF := i;
    END;
    rslt := JOIN(indy, dep, LEFT.rid=RIGHT.id AND LEFT.wi=RIGHT.wi, getD(LEFT,RIGHT),
                 LEFT OUTER, LIMIT(1,FAIL));
    RETURN rslt;
  END;  // to instance data
  // Model data
  SHARED Field_ID := 1;
  SHARED s_type_id:= Field_ID + 1;        // 2
  SHARED k_type_id:= s_type_id + 1;       // 3
  SHARED degree_id:= k_type_id + 1;       // 4
  SHARED gamma_id := degree_id + 1;       // 5
  SHARED coef0_id := gamma_id + 1;        // 6
  SHARED k_id := coef0_id + 1;            // 7
  SHARED l_id := k_id + 1;                // 8
  SHARED scale_id := l_id + 1;            // 9
  SHARED scaleInfo_id := scale_id + 1;    //10
  SHARED pairs_a_id := scaleInfo_id + 1;  //11
  SHARED pairs_b_id := pairs_a_id + 1;    //12
  SHARED n_label_id := pairs_b_id + 1;    //13
  SHARED n_sv_id := n_label_id + 1;       //14
  SHARED lf_id := n_sv_id; //last field id
  SHARED UNSIGNED4 get_SV_Count(DATASET(SVM.Types.SVM_SV) sv) := FUNCTION
    sv_count := COUNT(sv);
    features := SUM(sv, COUNT(sv.features));
    RETURN 2*sv_count + 2*features;
  END;
  SHARED SVM.Types.R8Entry norm_feature(SVM.Types.SVM_Feature f, UNSIGNED c):=TRANSFORM
    SELF.v := CHOOSE(c, (REAL8)f.nominal, f.v);
  END;
  SHARED Work1 := RECORD
    DATASET(SVM.Types.R8Entry) d;
  END;
  SHARED Work1 cvt_2_r8(SVM.Types.SVM_SV s) := TRANSFORM
    SELF.d := DATASET([{2*(COUNT(s.features)+1)},{s.v_ord}], SVM.Types.R8Entry)
            + NORMALIZE(s.features, 2, norm_feature(LEFT, COUNTER));
  END;
  SHARED SVM.Types.R8Entry cvt_sv(DATASET(SVM.Types.SVM_SV) sv) := FUNCTION
    w1 := PROJECT(sv, cvt_2_r8(LEFT));
    rslt := NORMALIZE(w1, LEFT.d, TRANSFORM(SVM.Types.R8Entry, SELF.v:=RIGHT.v));
    RETURN rslt;
  END;

  SHARED FromModelRow(UNSIGNED id_base = 1000, DATASET(SVM.Types.Model) mdl) := FUNCTION
    ML_Types.Layout_Model getField(SVM.Types.Model m, UNSIGNED c) := TRANSFORM
      SELF.wi := m.wi;
      SELF.id := id_base + m.id;
      SELF.number := c;
      SELF.value := CHOOSE(c, (REAL8) m.id, (REAL8) m.svmType,
                              (REAL8) m.kernelType, (REAL8) m.degree,
                              m.gamma, m.coef0, (REAL8) m.k, (REAL8) m.l,
                              (REAL8) m.scale, (REAL8) COUNT(m.scaleInfo)*2,
                              (REAL8) COUNT(m.probA), (REAL8) COUNT(m.probB),
                              (REAL8) COUNT(m.labels), (REAL8) get_SV_Count(m.sv));
    END;
    fixed_part := NORMALIZE(mdl, lf_id-field_id+1, getField(LEFT,COUNTER));
    ML_Types.Layout_Model normSV(SVM.Types.Model m, SVM.Types.R8Entry e, UNSIGNED8 c) := TRANSFORM
      SELF.wi := m.wi;
      SELF.id := id_base + m.id;
      SELF.number := lf_id + c;
      SELF.value  := e.v;
    END;
    sv_part := NORMALIZE(mdl, cvt_SV(LEFT.sv), normSV(LEFT, RIGHT, COUNTER));
    ML_Types.Layout_Model normMean(SVM.Types.Model m, SVM.Types.FeatureStats ft,
                                UNSIGNED8 c, UNSIGNED base) := TRANSFORM
      SELF.wi := m.wi;
      SELF.id := id_base + m.id;
      SELF.number := base + c;
      SELF.value  := ft.mean;
    END;
    ML_Types.Layout_Model normSD(SVM.Types.Model m, SVM.Types.FeatureStats ft,
                                UNSIGNED8 c, UNSIGNED base) := TRANSFORM
      SELF.wi := m.wi;
      SELF.id := id_base + m.id;
      SELF.number := base + c;
      SELF.value  := ft.sd;
    END;
    mean_part := NORMALIZE(mdl, LEFT.ScaleInfo,
                          normMean(LEFT, RIGHT, COUNTER,
                                 lf_id + get_SV_Count(LEFT.sv)));
    sd_part   := NORMALIZE(mdl, LEFT.ScaleInfo,
                          normSD(LEFT, RIGHT, COUNTER,
                                 lf_id + get_SV_Count(LEFT.sv) + COUNT(LEFT.ScaleInfo)));
    ML_Types.Layout_Model normR8(SVM.Types.Model m, SVM.Types.R8Entry d,
                                 UNSIGNED c, UNSIGNED4 base) := TRANSFORM
      SELF.wi := m.wi;
      SELF.id := id_base + m.id;
      SELF.number := base + c;
      SELF.value := d.v;
    END;
    coef_part := NORMALIZE(mdl, LEFT.sv_coef,
                          normR8(LEFT, RIGHT, COUNTER,
                                 lf_id + get_SV_Count(LEFT.sv)
                                 + 2*COUNT(LEFT.ScaleInfo)));
    rho_part  := NORMALIZE(mdl, LEFT.rho,
                          normR8(LEFT, RIGHT, COUNTER,
                                lf_id + get_SV_Count(LEFT.sv)
                                + 2*COUNT(LEFT.ScaleInfo)
                                + COUNT(LEFT.sv_coef)));
    probA_part:= NORMALIZE(mdl, LEFT.probA,
                           normR8(LEFT, RIGHT, COUNTER,
                                  lf_id + get_SV_Count(LEFT.sv)
                                  + 2*COUNT(LEFT.ScaleInfo)
                                  + COUNT(LEFT.sv_coef)
                                  + COUNT(LEFT.rho)));
    probB_part:= NORMALIZE(mdl, LEFT.probB,
                           normR8(LEFT, RIGHT, COUNTER,
                                  lf_id + get_SV_Count(LEFT.sv)
                                  + 2*COUNT(LEFT.ScaleInfo)
                                  + COUNT(LEFT.sv_coef)
                                  + COUNT(LEFT.rho)
                                  + COUNT(LEFT.probA)));
    ML_Types.Layout_Model normI4(SVM.Types.Model m, SVM.Types.I4Entry d,
                                 UNSIGNED c, UNSIGNED4 base) := TRANSFORM
      SELF.wi := m.wi;
      SELF.id := id_base + m.id;
      SELF.number := base + c;
      SELF.value := (REAL8) d.v;
    END;
    lab_part  := NORMALIZE(mdl, LEFT.labels,
                           normI4(LEFT, RIGHT, COUNTER,
                                  lf_id + get_SV_Count(LEFT.sv)
                                  + 2*COUNT(LEFT.ScaleInfo)
                                  + COUNT(LEFT.sv_coef)
                                  + COUNT(LEFT.rho)
                                  + COUNT(LEFT.probA)
                                  + COUNT(LEFT.probB)));
    nsv_part  := NORMALIZE(mdl, LEFT.nSV,
                           normI4(LEFT, RIGHT, COUNTER,
                                  lf_id + get_SV_Count(LEFT.sv)
                                  + 2*COUNT(LEFT.ScaleInfo)
                                  + COUNT(LEFT.sv_coef)
                                  + COUNT(LEFT.rho)
                                  + COUNT(LEFT.probA)
                                  + COUNT(LEFT.probB)
                                  + COUNT(LEFT.Labels)));
    RETURN fixed_part
         + sv_part
         + mean_part
         + sd_part
         + coef_part
         + rho_part
         + probA_part
         + probB_part
         + lab_part
         + nsv_part;
  END;

  /**
   * Convert from SVM Model type to standardized Layout_Model format. The
   * Layout_Model format is harder to interpret, but more generalized.
   * @param id_base Base number from which to start model IDs (default: 1000).
   * @param mdl Object of SVM Model type (see SupportVectorMachines.Types).
   * @return Convert SVM model in Layout_Model format (see ML_Core.Types for
   * format definition).
   */
  EXPORT FromModel(UNSIGNED id_base = 1000, DATASET(SVM.Types.Model) mdl) := FUNCTION
    {UNSIGNED2 wi, INTEGER4 id, DATASET(ML_Types.Layout_Model) mdls} convertModel(SVM.Types.Model mdl) := TRANSFORM
      SELF.wi := mdl.wi;
      SELF.id := mdl.id;
      SELF.mdls := FromModelRow(id_base, PROJECT(DATASET(mdl), SVM.Types.Model));
    END;
    mdls_LM := PROJECT(mdl, convertModel(LEFT));
    mdls_LM_sort := SORT(mdls_LM, wi, id);

    ML_Types.Layout_Model combineMdls(ML_Types.Layout_Model mdl_LM) := TRANSFORM
      SELF := mdl_LM;
    END;
    rslt := NORMALIZE(mdls_LM, LEFT.mdls, combineMdls(RIGHT));
    RETURN rslt;
  END;



  SHARED Exp_NF := RECORD(ML_Types.NumericField)
    UNSIGNED feature;
    UNSIGNED vector;
  END;
  SHARED Exp_NF mark_xnf(Exp_NF prev, Exp_NF curr) := TRANSFORM
    SELF.vector := IF(curr.number>prev.vector,
                      curr.number+(UNSIGNED)curr.value-1,
                      prev.vector);
    SELF.feature:= IF(curr.number=prev.feature, prev.feature, curr.number+1);
    SELF := curr;
  END;
  SHARED SVM.Types.SVM_Feature makeFeature(DATASET(Exp_NF) d) := TRANSFORM
    SELF.nominal := (INTEGER) d[1].value;
    SELF.v := d[2].value;
  END;
  SHARED SVM.Types.SVM_SV roll_xnf(Exp_NF f, DATASET(Exp_NF) ds) := TRANSFORM
    feature_data := GROUP(CHOOSEN(ds, ALL, 3), feature);
    SELF.v_ord := (INTEGER) ds[2].value;
    SELF.features := ROLLUP(feature_data, GROUP, makeFeature(ROWS(LEFT)));
  END;
  SHARED SVM.Types.R8Entry toR8(ML_Types.Layout_Model nf) := TRANSFORM
    SELF.v := nf.value;
  END;
  SHARED SVM.Types.I4Entry toI4(ML_Types.Layout_Model nf) := TRANSFORM
    SELF.v := (INTEGER) nf.value;
  END;
  SHARED toSV_Dataset(DATASET(ML_Types.Layout_Model) nf) := FUNCTION
    x_nf := PROJECT(nf, TRANSFORM(Exp_NF, SELF:=LEFT,SELF:=[]));
    marked_x := ITERATE(x_nf, mark_xnf(LEFT,RIGHT));
    group_x := GROUP(marked_x, vector);
    rslt := ROLLUP(group_x, GROUP, roll_xnf(LEFT, ROWS(LEFT)));
    RETURN rslt;
  END;
  SHARED toFeatStat_Dataset(DATASET(ML_Types.Layout_Model) nf) := FUNCTION
    statCount := COUNT(nf)/2;
    numberStart := MIN(nf, number);
    means := nf(number BETWEEN numberStart AND numberStart+statCount);
    sds := nf(number BETWEEN numberStart+statCount+1 AND numberStart+statCount*2);

    SVM.Types.FeatureStats makeStats(ML_Types.Layout_Model m, ML_Types.Layout_Model sd)
    := TRANSFORM
      SELF.indx := IF(
        m.number = numberStart, -1,
        m.number - numberStart + 1);
      SELF.mean := m.value;
      SELF.sd := sd.value;
    END;

    rslt := JOIN(means, sds, LEFT.number = RIGHT.number-statCount, makeStats(LEFT, RIGHT));
    RETURN rslt;
  END;

  /**
   * Convert from standardized Layout_Model format (see ML_Core.Types) to
   * SVM Model format (see SupportVectorMachines.Types). The SVM model format
   * is less general, but easier to interpret.
   * @param mdl Trained SVM in Layout_Model format.
   * @return Converted SVM model in SVM Model format.
   */
  EXPORT ToModel(DATASET(ML_Types.Layout_Model) mdl) := FUNCTION
    mdl_grp := GROUP(SORTED(mdl, wi, id, number, ASSERT), wi, id);

    SVM.Types.Model rollModel(DATASET(ML_Types.Layout_Model) d) := TRANSFORM
      fixed := DICTIONARY(d(number<=lf_id), {number=>value});
      nsv := (UNSIGNED) fixed[n_sv_id].value;
      probA := (UNSIGNED) fixed[pairs_a_id].value;
      probB := (UNSIGNED) fixed[pairs_b_id].value;
      labels:= (UNSIGNED) fixed[n_label_id].value;
      k := (UNSIGNED) fixed[k_id].value;
      l := (UNSIGNED) fixed[l_id].value;
      scale := (BOOLEAN) fixed[scale_id].value;
      scaleInfo := (UNSIGNED) fixed[scaleInfo_id].value;
      sv_start := lf_id + 1;
      sv_stop := sv_start + nsv - 1;
      scaleInfo_start := sv_stop + 1;
      scaleInfo_stop := scaleInfo_start + scaleInfo - 1;
      coef_start := scaleInfo_stop + 1;
      coef_stop := coef_start + (k-1)*l - 1;
      rho_start := coef_stop + 1;
      rho_stop := rho_start + (k-1)*k/2 - 1;
      probA_start := rho_stop + 1;
      probA_stop := probA_start + probA - 1;
      probB_start := probA_stop + 1;
      probB_stop := probB_start + probB - 1;
      label_start := probB_stop + 1;
      label_stop := label_start + labels - 1;
      nSV_start := label_stop + 1;
      nSV_stop := nSV_start + labels - 1;
      SELF.wi := d[1].wi;
      SELF.id := (UNSIGNED) fixed[Field_ID].value;
      SELF.svmType := (UNSIGNED1) fixed[s_type_id].value;
      SELF.kernelType := (UNSIGNED1) fixed[k_type_id].value;
      SELF.degree := (UNSIGNED) fixed[degree_id].value;
      SELF.coef0 := fixed[coef0_id].value;
      SELF.gamma := fixed[gamma_id].value;
      SELF.k := k;
      SELF.l := l;
      SELF.scale := scale;
      SELF.scaleInfo := toFeatStat_Dataset(d(number BETWEEN scaleInfo_start AND scaleInfo_stop));
      SELF.sv := toSV_Dataset(d(number BETWEEN sv_start AND sv_stop));
      SELF.sv_coef := PROJECT(d(number BETWEEN coef_start AND coef_stop), toR8(LEFT));
      SELF.rho := PROJECT(d(number BETWEEN rho_start AND rho_stop), toR8(LEFT));
      SELF.probA := PROJECT(d(number BETWEEN probA_start AND probA_stop), toR8(LEFT));
      SELF.probB := PROJECT(d(number BETWEEN probB_start AND probB_stop), toR8(LEFT));
      SELF.labels:= PROJECT(d(number BETWEEN label_start AND label_stop), toI4(LEFT));
      SELF.nSV := PROJECT(d(number BETWEEN nSV_start AND nSV_stop), toI4(LEFT));
    END;
    rslt := ROLLUP(SORTED(mdl_grp, wi, id, number, ASSERT), GROUP, rollModel(ROWS(LEFT)));
    RETURN rslt;
  END;

  // Model := SVM.Types.Model;
  // Model_Split := RECORD
    // Model.wi;
    // Model.id;
    // Model.svmType;
    // Model.kernelType;
    // Model.degree;
    // Model.gamma;
    // Model.coef0;
    // Model.k;
    // Model.l;
    // Model.scale;
    // UNSIGNED8 scaleInfo_cnt;
    // UNSIGNED8 probA_cnt;
    // UNSIGNED8 probB_cnt;
    // UNSIGNED8 label_cnt;
    // UNSIGNED8 sv_cnt;
    // DATASET(ML_Types.Layout_Model) sv_LM;
    // DATASET(ML_Types.Layout_Model) scaleInfo_LM;
    // DATASET(ML_Types.Layout_Model) sv_coef_LM;
    // DATASET(ML_Types.Layout_Model) rho_LM;
    // DATASET(ML_Types.Layout_Model) probA_LM;
    // DATASET(ML_Types.Layout_Model) probB_LM;
    // DATASET(ML_Types.Layout_Model) labels_LM;
    // DATASET(ML_Types.Layout_Model) nSV_LM;
  // END;

  // EXPORT ToModel(DATASET(ML_Types.Layout_Model) mdl) := FUNCTION
    // Model_Split getParentRows(ML_Types.Layout_Model mdl) := TRANSFORM
      // SELF := mdl;
      // SELF := [];
    // END;
    // mdlParent := PROJECT(mdl, getParentRows(LEFT));
    // Model_Split denormLM(Model_Split p, ML_Types.Layout_Model c) := TRANSFORM
      // SELF.wi            := c.wi;
      // SELF.id            := IF(c.number = Field_ID, (INTEGER4) c.value, p.id);
      // SELF.svmType       := IF(c.number = s_type_id, (UNSIGNED1) c.value, p.svmType);
      // SELF.kernelType    := IF(c.number = k_type_id, (UNSIGNED1) c.value, p.kernelType);
      // SELF.degree        := IF(c.number = degree_id, (UNSIGNED) c.value, p.degree);
      // SELF.coef0         := IF(c.number = coef0_id, c.value, p.coef0);
      // SELF.gamma         := IF(c.number = gamma_id, c.value, p.gamma);
      // SELF.k             := IF(c.number = k_id, (UNSIGNED) c.value, p.k);
      // SELF.l             := IF(c.number = l_id, (UNSIGNED) c.value, p.l);
      // SELF.scale         := IF(c.number = scale_id, (BOOLEAN) c.value, p.scale);
      // SELF.scaleInfo_cnt := IF(c.number = scaleInfo_id, (UNSIGNED) c.value, p.scaleInfo_cnt);
      // SELF.probA_cnt     := IF(c.number = pairs_a_id, (UNSIGNED) c.value, p.probA_cnt);
      // SELF.probB_cnt     := IF(c.number = pairs_b_id, (UNSIGNED) c.value, p.probB_cnt);
      // SELF.label_cnt     := IF(c.number = n_label_id, (UNSIGNED) c.value, p.label_cnt);
      // SELF.sv_cnt        := IF(c.number = n_sv_id, (UNSIGNED) c.value, p.sv_cnt);
      // SELF.sv_LM         := IF(c.number BETWEEN lf_id + 1 AND lf_id + p.sv_cnt,
        // p.sv_LM + c, p.sv_LM);
        // SELF.scaleInfo_LM := IF(c.number BETWEEN
        // lf_id + p.sv_cnt + 1 AND
        // lf_id + p.sv_cnt + p.scaleInfo_cnt,
        // p.scaleInfo_LM + c,
        // p.scaleInfo_LM);
        // SELF.sv_coef_LM := IF(c.number BETWEEN
        // lf_id + p.sv_cnt + p.scaleInfo_cnt + 1 AND
        // lf_id + p.sv_cnt + p.scaleInfo_cnt + (p.k-1)*p.l,
        // p.sv_coef_LM + c,
        // p.sv_coef_LM);
        // SELF.rho_LM := IF(c.number BETWEEN
        // lf_id + p.sv_cnt + p.scaleInfo_cnt + (p.k-1)*p.l + 1 AND
        // lf_id + p.sv_cnt + p.scaleInfo_cnt + (p.k-1)*(p.l + p.k/2),
        // p.rho_LM + c,
        // p.rho_LM);
        // SELF.probA_LM := IF(c.number BETWEEN
        // lf_id + p.sv_cnt + p.scaleInfo_cnt + (p.k-1)*(p.l + p.k/2) + 1 AND
        // lf_id + p.sv_cnt + p.scaleInfo_cnt + (p.k-1)*(p.l + p.k/2) + p.probA_cnt,
        // p.probA_LM + c,
        // p.probA_LM);
        // SELF.probB_LM := IF(c.number BETWEEN
        // lf_id + p.sv_cnt + p.scaleInfo_cnt + (p.k-1)*(p.l + p.k/2) + p.probA_cnt + 1 AND
        // lf_id + p.sv_cnt + p.scaleInfo_cnt + (p.k-1)*(p.l + p.k/2) + p.probA_cnt + p.probB_cnt,
        // p.probB_LM + c,
        // p.probB_LM);
        // SELF.labels_LM := IF(c.number BETWEEN
        // lf_id + p.sv_cnt + p.scaleInfo_cnt + (p.k-1)*(p.l + p.k/2) + p.probA_cnt + p.probB_cnt + 1 AND
        // lf_id + p.sv_cnt + p.scaleInfo_cnt + (p.k-1)*(p.l + p.k/2) + p.probA_cnt + p.probB_cnt + p.label_cnt,
        // p.labels_LM + c,
        // p.labels_LM);
        // SELF.nSV_LM := IF(c.number BETWEEN
        // lf_id + p.sv_cnt + p.scaleInfo_cnt + (p.k-1)*(p.l + p.k/2) + p.probA_cnt + p.probB_cnt + p.label_cnt + 1 AND
        // lf_id + p.sv_cnt + p.scaleInfo_cnt + (p.k-1)*(p.l + p.k/2) + p.probA_cnt + p.probB_cnt + p.label_cnt + p.label_cnt,
        // p.nSV_LM + c,
        // p.nSV_LM);
      // SELF := p;
    // END;
    // mdlLumped := DENORMALIZE(
      // mdlParent,
      // mdl,
      // LEFT.wi = RIGHT.wi AND LEFT.id = RIGHT.id,
      // denormLM(LEFT, RIGHT),
      // NOSORT
    // );
    // Model getSVMModel(Model_Split mL) := TRANSFORM
      // SELF.scaleInfo := toFeatStat_Dataset(mL.scaleInfo_LM);
      // SELF.sv := toSV_Dataset(mL.sv_LM);
      // SELF.sv_coef := PROJECT(mL.sv_coef_LM, toR8(LEFT));
      // SELF.rho := PROJECT(mL.rho_LM, toR8(LEFT));
      // SELF.probA := PROJECT(mL.probA_LM, toR8(LEFT));
      // SELF.probB := PROJECT(mL.probB_LM, toR8(LEFT));
      // SELF.labels := PROJECT(mL.labels_LM, toI4(LEFT));
      // SELF.nSV := PROJECT(mL.nSV_LM, toI4(LEFT));
      // SELF := mL;
    // END;
    // rslt := PROJECT(mdlLumped, getSVMModel(LEFT));
    // RETURN rslt;
  // END;
END;