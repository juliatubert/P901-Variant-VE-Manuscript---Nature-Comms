/*********************************************************************************************************************************************************
Program: Analysis_2nd_variant_VE

Purpose: 
    Moderna 2nd variant vaccine effectiveness (VE) analysis: 
    Describe the study population, baseline characterizes 
    and conduct the analysis for study outcomes.

Note: Table and figure #s may differ in actual publication 
*********************************************************************************************************************************************************/

*** Set up infile passway;
libname in "~/2nd_variant_VE" access=readonly;
%let infile = in.analysis_2nd_variant_VE_revised;

*** Set up output passway;
%let outpath = 2nd_variant_VE;
libname out  "~/&outpath";
%let gpath = "~/&outpath"; 
%let update_date=%sysfunc(today(), date9.);

*** Call macros;
%include "~/&outpath./Analysis_2nd_variant_VE_macro.sas";

* Prepare analysis datasets;
    data analysis_data; 
	length days_interval variant variant_var variant_sgtf variant_var_sgtf $15;
	set &infile;
	* Keep matched only;
	where matched=1;
	* Edit days_interval;
	if days_interval='         .' then days_interval='';
	if 91<=time_interval<=150 then days_interval='91-150 days';
	if time_interval>150 then days_interval='>150 days';
	* Create new interval for time between 2nd/3rd dose and test date;
	if .<vac_date2<index_date then time_interval_dose2=index_date-vac_date2;
	if .<vac_date3<index_date then time_interval_dose3=index_date-vac_date3;
	* Create vac_status variable;
	if arm='Unvaccinated' then vac_status='Unvaccinated';
	else vac_status=cohort;
    * Create vac_status_var variable (used for var names, can't start with number);
    if arm='Unvaccinated'   then vac_status_var='Unvaccinated';
    else if cohort='1-dose' then vac_status_var='Vac_1dose'; 
    else if cohort='2-dose' then vac_status_var='Vac_2dose';
    else if cohort='3-dose' then vac_status_var='Vac_3dose'; 
    else if cohort='4-dose' then vac_status_var='Vac_4dose';
    * Create variant var (BA.2* excludes BA.2.12.1);
    if variant_label='BA.2 except' then variant='BA.2*'; 
	else variant=variant_label;
	* Create variant vars with SGTF data for sensitivity analysis;
	if results_sgtf^='' then sgtf_available=1; else sgtf_available=0;
	variant_sgtf = variant;
	if variant in ('BA.4','BA.5') then variant_SGTF='BA.4/BA.5';
	if variant in ('BA.2*','BA.2.12.1') then variant_SGTF='BA.2';
	if qc_status='fail' then do;
		if index_month in ('2022-01', '2022-02', '2022-03', '2022-04') and results_sgtf='SGTF+' then variant_sgtf='BA.1';
		if index_month in ('2022-05', '2022-06') and results_sgtf='SGTF+' then variant_sgtf='BA.4/BA.5';
		if index_month='2022-01' and results_sgtf='SGTF-' then variant_sgtf='Delta';
		if index_month^='2022-01' and results_sgtf='SGTF-' then variant_sgtf='BA.2';
		end;
	* Create variant_var (used for var names);
	variant_var = compress(variant, '/.,* ');
	variant_var_sgtf = compress(variant_sgtf, '/.,* ');
	run; 

**********************************
** Table 1. Sequencing characteristics of SARS-CoV-2 specimens, by sequencing status and mRNA-1273 vaccination;

* Get frequencies and percentages;
 	data table1data; set analysis_data;
    where COVID19='P';
	run; 
    data success fail; set table1data;
    if qc_status='pass' then output success;
    else output fail;
    run;
    * For sequencing success part;
    %stat_TABLE(DSN=success, ID=mrn, by=vac_status_var,
                var=Specimen_Type ct_value sgtf_available covid_hosp covid_hosp_death variant, 
                type=2 2 2 , pvalues=N, 
                PCTPRSNT=2, outdat=Table1Part1);
    * For sequencing failure part;
    %stat_TABLE(DSN=fail, ID=mrn, by=vac_status_var,
                var=Specimen_Type ct_value sgtf_available covid_hosp covid_hosp_death, 
                type=2 2, pvalues=N, 
                PCTPRSNT=2, outdat=Table1Part2);

* Minor adjustment of the macro output to meet SAP table shell structure and format;
    data Table1Part1; set Table1Part1;
    order=_N_;
    * Create a total line at begining;
    if _N_=1 then do; level='Total'; 
	c1=compress(put(byn1,$10.),' '); 
	c2=compress(put(byn2,$10.),' '); 
	c3=compress(put(byn3,$10.),' '); 
	c4=compress(put(byn4,$10.),' '); 
	c5=compress(put(byn5,$10.),' '); 
	ctotal=compress(put(nobs,$10.),' '); 
    end;
    run;
    data Table1Part2; set Table1Part2;
    order=_N_;
    * Create a total line at begining;
    if _N_=1 then do; level='Total'; 
	c1=compress(put(byn1,$10.),' '); 
	c2=compress(put(byn2,$10.),' '); 
	c3=compress(put(byn3,$10.),' '); 
	c4=compress(put(byn4,$10.),' '); 
	c5=compress(put(byn5,$10.),' '); 
	ctotal=compress(put(nobs,$10.),' '); 
    end;
    run;
	* Modify order<# based on row number in Table1part1;
    proc sql; create table Table1 as
    select a.level, 
           a.c2 as succ_1dose, a.c3 as succ_2dose, a.c4 as succ_3dose, a.c5 as succ_4dose, a.c1 as succ_unv, a.ctotal as succ_total, ' ' as separator1,
           case when a.order<24 then b.c2 else 'N/A' end as fail_1dose, 
           case when a.order<24 then b.c3 else 'N/A' end as fail_2dose, 
           case when a.order<24 then b.c4 else 'N/A' end as fail_3dose, 
           case when a.order<24 then b.c5 else 'N/A' end as fail_4dose, 
           case when a.order<24 then b.c1 else 'N/A' end as fail_unv, 
           case when a.order<24 then b.ctotal else 'N/A' end as fail_total
      from Table1part1 a left join Table1part2 b on a.order=b.order 
     where a.level not in ('','sgtf_available','covid_hosp','covid_hosp_death','    0') /*only keep level=1 for binary vars*/
     order by a.order;
    quit;

	* Confirm total counts;
	proc freq data=table1data;
		table qc_status*vac_status/norow nocol nopercent;
		table sgtf_available*qc_status;
	run;

	title 'Table 1';
    proc print data=Table1; run;
	title;

**********************************
** Figure 4. Distribution of SARS-CoV-2 variants by mRNA-1273 vaccination status;

proc freq data=table1data order=data noprint;
	table vac_status*variant / sparse list out=FreqOutSortedf4;
	where variant^='Unidentified';
run;
proc freq data=table1data noprint;
    tables vac_status / out=FreqOutlabelf4(drop=percent);
	where variant^='Unidentified';
run;
proc sql; create table figure4data as
	select cat(a.vac_status,repeat(' ',15-int(log10(b.count))),'(N=', b.count,')') as Xlabel, a.count, (100*a.count/sum(a.count)) as percent, a.variant
    from FreqOutSortedf4 a left join FreqOutlabelf4 b on a.vac_status=b.vac_status
	where variant^='Unidentified'
	group by a.vac_status
    order by a.vac_status, a.variant;
quit;

* Generate plot;
title 'Figure 4';
ods graphics /  noborder imagename="Figure4_&update_date";
ods listing gpath=&gpath. image_dpi=300;
proc sgplot data=figure4data;
	styleattrs datacolors=( cx619070 /*BA.1 - dusty green*/
						    cxe5b82e /*BA.2 - yellow*/ 
							cxff8000 /*BA.2.12.1 - orange*/
							cx82d0c9 /*BA.3 - light turquoise*/
							cx99293d /*BA.4 - dark pink*/
						    cxd47676 /*BA.5 - burnt salmon*/
					    	cx8b6e8b /*Delta - dusty purple*/
							cx143332 /*Recombinant - deep bluish green*/
						    );
	vbar Xlabel / response=Percent group=variant 
	              grouporder=ascending groupdisplay=stack nooutline; 
	xaxis discreteorder=data label='Vaccination status' labelattrs=(size = 12); 
	yaxis grid values=(0 to 100 by 10) label="Percentage" labelattrs=(size = 12);
	keylegend / title="";
run;
title;


**********************************
** Figure 5. Distribution of SARS-CoV-2 variants by month of sample collection;

proc freq data=table1data order=data noprint;
	table index_date*variant / sparse list out=FreqOutSortedf5;
	format index_date monyy7.;
	where variant^='Unidentified';
run;
proc freq data=table1data noprint;
    tables index_date / out=FreqOutlabelf5(drop=percent);
	format index_date monyy7.;
	where variant^='Unidentified';
run;
proc sql; create table figure5data as
    select cat(put(a.index_date,monyy7.),repeat(' ',15-int(log10(b.count))),'(N=', b.count,')') as Xlabel, a.count, (100*a.count/sum(a.count)) as percent, a.variant
    from FreqOutSortedf5 a left join FreqOutlabelf5 b on a.index_date=b.index_date
	where variant^='Unidentified'
	group by a.index_date
    order by a.index_date, a.variant;
quit;

* Generate plot;
title 'Figure 5';
ods graphics /  noborder imagename="Figure5_&update_date";
ods listing gpath=&gpath. image_dpi=300;
proc sgplot data=figure5data;
	styleattrs datacolors=( cx619070 /*BA.1 - dusty green*/
						    cxe5b82e /*BA.2 - yellow*/ 
							cxff8000 /*BA.2.12.1 - orange*/
							cx82d0c9 /*BA.3 - light turquoise*/
							cx99293d /*BA.4 - dark pink*/
						    cxd47676 /*BA.5 - burnt salmon*/
					    	cx8b6e8b /*Delta - dusty purple*/
							cx143332 /*Recombinant - deep bluish green*/
							);
	vbar Xlabel / response=Percent group=variant 
	              grouporder=ascending groupdisplay=stack nooutline; 
	xaxis discreteorder=data label='Month of specimen collection' labelattrs=(size = 12); 
	yaxis grid values=(0 to 100 by 10) label="Percentage" labelattrs=(size = 12);
	keylegend / title="";
run;
title;


**********************************
** Table 3. Characteristics of SARS-CoV-2 test-positive cases and test-negative controls;

%TNDbaselinetable(analysis_data);
title 'Table 3';
proc print data=Table3;
title;

* View covariates with p<0.01 and ASD>0.1 for model adjustment;
    data significant; set Table3;
    pvalue=.; pvalue=p; ABS_value=.; ABS_value=ABS;
    if (.<pvalue<0.1 or p='<0.01') and ABS>0.1 and var not in ('frailty_index','age','vac_status'); 
    order=_N_;
    keep var p ABS order;
    run;
	* hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month (already adjusting for matching vars: age_grp sex Race_eth index_month); 


**********************************
** Table 5. Odds ratio of test-positive SARS-CoV-2 cases;

* Create subgroupflag for variant subgroup (matching sets);
    proc sql; create table subgroupflag as
    select a.*, b.variant_var as variant_subgroup, b.variant_var_sgtf as variant_subgroup_sgtf
    from analysis_data a left join table1data b on a.pair_id=b.pair_id
	order by pair_id, covid19;
    quit; 
	* Create collapsed versions of covariates for convergence issues;
	data subgroupflag; set subgroupflag;
	* Index_month;
	if index_month in ('2022-03','2022-04') then index_month_MarApr='Mar-Apr'; else index_month_MarApr=index_month;
	if index_month in ('2022-04','2022-05') then index_month_AprMay='Apr-May'; else index_month_AprMay=index_month;
	if index_month in ('2022-03','2022-04','2022-05') then index_month_MarAprMay='Mar-May'; else index_month_MarAprMay=index_month;
	* VA_grp;
	length VA_grp_3grps VA_grp_2grps $5;
	if VA_grp in ('0','1-4') then VA_grp_3grps='0-4'; else VA_grp_3grps=VA_grp;
	if VA_grp in ('0','1-4','5-10') then VA_grp_2grps='0-10'; else VA_grp_2grps=VA_grp;
	run;
    *QC;
    proc freq data=subgroupflag;
	table variant_subgroup*vac_status*COVID19/list missing;
    table VA_grp*VA_grp_3grps*VA_grp_2grps/list missing;
	table index_month*index_month_MarApr*index_month_AprMay*index_month_MarAprMay/list missing;
    run;

/*** 3-dose vs. unvac: VE against covid_dx ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=3-dose;
	%let comparison_grp=Unvaccinated;
	%let suffix=3dose_unvac;
	%let vars_num=;
	* BA.1; 
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA1, 	condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA1, 	condition_model=0, time_since_vac=1); 
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2, 	condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2,	condition_model=0, time_since_vac=1);
	* BA.2.12.1;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2121, condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2121, condition_model=0, time_since_vac=1);
	* BA.4;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month;
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4, 	condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4, 	condition_model=0, time_since_vac=1);
	* BA.5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA5, 	condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA5, 	condition_model=0, time_since_vac=1);
	* Combine results;
    %TND_combine(table5_&suffix);

/*** 3-dose vs. 2-dose: VE against covid_dx  ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=3-dose;
	%let comparison_grp=2-dose;
	%let suffix=3dose_2dose;
	%let vars_num=time_interval_dose2;
	* BA.1; 
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA1, 	condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA1, 	condition_model=0, time_since_vac=1); 
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2, 	condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2,	condition_model=0, time_since_vac=1);
	* BA.2.12.1;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2121, condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2121, condition_model=0, time_since_vac=1);
	* BA.4;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month;
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4, 	condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4, 	condition_model=0, time_since_vac=1);
	* BA.5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month;  
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA5, 	condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA5, 	condition_model=0, time_since_vac=1);
	* Combine results;
    %TND_combine(table5_&suffix);

/*** 4-dose vs. 3-dose: VE against covid_dx  ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=4-dose;
	%let comparison_grp=3-dose;
	%let suffix=4dose_3dose;
	%let vars_num=time_interval_dose3;
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2, 	condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2,	condition_model=0, time_since_vac=1);
	* BA.2.12.1;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2121, condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2121, condition_model=0, time_since_vac=1);
	* BA.4;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4, 	condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4, 	condition_model=0, time_since_vac=1);
	* BA.5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA5, 	condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA5, 	condition_model=0, time_since_vac=1);
	* Combine results;
    %TND_combine(table5_&suffix);

/*** 4-dose vs. unvac: VE against covid_dx  ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=4-dose;
	%let comparison_grp=Unvaccinated;
	%let suffix=4dose_unvac;
	%let vars_num=;
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2, 	condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2,	condition_model=0, time_since_vac=1);
	* BA.2.12.1;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month_MarApr; *Collapsed index_month;
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2121, condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2121, condition_model=0, time_since_vac=1);
	* BA.4;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month_AprMay; *Collapsed index_month;
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4, 	condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4, 	condition_model=0, time_since_vac=1);
	* BA.5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA5, 	condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA5, 	condition_model=0, time_since_vac=1);
	* Combine results;
    %TND_combine(table5_&suffix);


**********************************
** Table 6. Odds ratio of test-positive hosptializated SARS-CoV-2 cases;

* Remove COVID infection sets, the infection is not considered case or control here;
    proc sql; create table hospitalsets_tmp as
    select distinct a.* from subgroupflag a, subgroupflag(where=(covid_hosp)) b
     where a.pair_id=b.pair_id;
    quit;

* Check # of outcome available;
    proc freq data=hospitalsets_tmp;
    where covid_hosp;
    table vac_status*variant_subgroup /list missing;
    run;

* Combine variant_subgroup for hosp analyses, due to small counts;
    data hospitalsets; set hospitalsets_tmp;
	if variant_subgroup in ('BA2','BA2121') then variant_subgroup='BA2';
	if variant_subgroup in ('BA4','BA5') then variant_subgroup='BA4BA5';
    run;
    proc freq data=hospitalsets;
	where covid_hosp;
    table vac_status*variant_subgroup/list missing;
	run;
    proc freq data=hospitalsets;
	where covid_hosp and variant_subgroup^='Unidentified' and vac_status in ('3-dose','4-dose');
    table vac_status*variant_subgroup*days_interval/list missing;
	run;
	proc sort data=hospitalsets; by vac_status variant_subgroup time_interval; run;
	proc print data=hospitalsets;
	var vac_status variant_subgroup time_interval;
	where covid_hosp and variant_subgroup^='Unidentified' and vac_status in ('3-dose','4-dose');
	run;


/*** 3-dose vs. unvac: VE against covid_hosp ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=3-dose;
	%let comparison_grp=Unvaccinated;
	%let suffix=3dose_unvac;
	%let vars_num=;
	* BA.1; 
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month_MarAprMay; *Collapsed index_month;
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA1, condition_model=0); 
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month; 
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA2, condition_model=0);
	* BA.4/5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month; 
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA4BA5, condition_model=0);
	* Combine results;
    %TND_combine(table6_&suffix);

/*** 3-dose vs. 2-dose: VE against covid_hosp ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=3-dose;
	%let comparison_grp=2-dose;
	%let suffix=3dose_2dose;
	%let vars_num=time_interval_dose2;
	* BA.1; 
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month_MarAprMay; *Collapsed index_month;
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA1, condition_model=0); 
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month; 
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA2, condition_model=0);
	* BA.4/5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month; 
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA4BA5, condition_model=0);
	* Combine results;
    %TND_combine(table6_&suffix);

/*** 4-dose vs. 3-dose: VE against covid_hosp ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=4-dose;
	%let comparison_grp=3-dose;
	%let suffix=4dose_3dose;
	%let vars_num=time_interval_dose3;
	* BA.2;
	%let vars_class=/*hist_covid_diag*/ hist_covid_test VA_grp_2grps /*medical_center_and_area*/ index_month; *Collapsed VA_grp further;
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA2, condition_model=0); 
	* BA.4/5;
	%let vars_class=/*hist_covid_diag*/ hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month; 
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA4BA5, condition_model=0); 
	* Combine results;
    %TND_combine(table6_&suffix);

/*** 4-dose vs. unvac: VE against covid_hosp ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=4-dose;
	%let comparison_grp=Unvaccinated;
	%let suffix=4dose_unvac;
	%let vars_num=;
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month_MarApr; *Collapsed index_month;
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA2, condition_model=0); 
	* BA.4/5;
	data hospitalsets2; set hospitalsets;
	if race_eth in ('Non-Hispanic Asian','Non-Hispanic Black','Other/Unknown') then race_eth='Other/Unknown';
	run;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month; 
	%TND(hospitalsets2,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA4BA5, condition_model=0); *Collapsed race_eth in hospitalsets2;
	* Combine results;
    %TND_combine(table6_&suffix);


**********************************
** Table 8-9. Sensitivity analysis with imputed variants based on SGTF data;

** SGTF sensitivity analyses for covid_dx;

/*** 3-dose vs. unvac: VE against covid_dx ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=3-dose;
	%let comparison_grp=Unvaccinated;
	%let suffix=3dose_unvac;
	%let vars_num=;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month; 
	* BA.1; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA1, variant_var=variant_subgroup_sgtf, condition_model=0); 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA1, variant_var=variant_subgroup_sgtf, condition_model=0, time_since_vac=1); 
	* BA.2;
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2, variant_var=variant_subgroup_sgtf, condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2, variant_var=variant_subgroup_sgtf, condition_model=0, time_since_vac=1); 
	* BA.4/5;
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4BA5, variant_var=variant_subgroup_sgtf, condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4BA5, variant_var=variant_subgroup_sgtf, condition_model=0, time_since_vac=1); 
	* Combine results;
    %TND_combine(table8_&suffix);


/*** 3-dose vs. 2-dose: VE against covid_dx ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=3-dose;
	%let comparison_grp=2-dose;
	%let suffix=3dose_2dose;
	%let vars_num=time_interval_dose2;
	* BA.1; 
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA1, variant_var=variant_subgroup_sgtf, condition_model=0); 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA1, variant_var=variant_subgroup_sgtf, condition_model=0, time_since_vac=1); 
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2, variant_var=variant_subgroup_sgtf, condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2, variant_var=variant_subgroup_sgtf, condition_model=0, time_since_vac=1); 
	* BA.4/5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4BA5, variant_var=variant_subgroup_sgtf, condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4BA5, variant_var=variant_subgroup_sgtf, condition_model=0, time_since_vac=1); 
	* Combine results;
    %TND_combine(table8_&suffix);

/*** 4-dose vs. 3-dose: VE against covid_dx ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=4-dose;
	%let comparison_grp=3-dose;
	%let suffix=4dose_3dose;
	%let vars_num=time_interval_dose3;
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2, variant_var=variant_subgroup_sgtf, condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2, variant_var=variant_subgroup_sgtf, condition_model=0, time_since_vac=1); 
	* BA.4/5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4BA5, variant_var=variant_subgroup_sgtf, condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4BA5, variant_var=variant_subgroup_sgtf, condition_model=0, time_since_vac=1); 
	* Combine results;
    %TND_combine(table8_&suffix);

/*** 4-dose vs. unvac: VE against covid_dx ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=4-dose;
	%let comparison_grp=Unvaccinated;
	%let suffix=4dose_unvac;
	%let vars_num=;
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2, variant_var=variant_subgroup_sgtf, condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2, variant_var=variant_subgroup_sgtf, condition_model=0, time_since_vac=1); 
	* BA.4/5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month_AprMay; *Collapsed index_month;
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4BA5, variant_var=variant_subgroup_sgtf, condition_model=0);
	%TND(subgroupflag,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4BA5, variant_var=variant_subgroup_sgtf, condition_model=0, time_since_vac=1); 
	* Combine results;
    %TND_combine(table8_&suffix);


** SGTF sensitivity analyses for covid_hosp;

/*** 3-dose vs. unvac: VE against covid_hosp ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=3-dose;
	%let comparison_grp=Unvaccinated;
	%let suffix=3dose_unvac;
	%let vars_num=;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month; 
	* BA.1; 
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA1, variant_var=variant_subgroup_sgtf, condition_model=0); 
	* BA.2;
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA2, variant_var=variant_subgroup_sgtf, condition_model=0);
	* BA.4/5;
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA4BA5, variant_var=variant_subgroup_sgtf, condition_model=0);
	* Combine results;
    %TND_combine(table9_&suffix);

/*** 3-dose vs. 2-dose: VE against covid_hosp ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=3-dose;
	%let comparison_grp=2-dose;
	%let suffix=3dose_2dose;
	%let vars_num=time_interval_dose2;
	* BA.1; 
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month_MarAprMay; *Collapsed index_month; 
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA1, variant_var=variant_subgroup_sgtf, condition_model=0); 	
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month; 
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA2, variant_var=variant_subgroup_sgtf, condition_model=0);
	* BA.4/5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month; 
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA4BA5, variant_var=variant_subgroup_sgtf, condition_model=0);
	* Combine results;
    %TND_combine(table9_&suffix);

/*** 4-dose vs. 3-dose: VE against covid_hosp ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=4-dose;
	%let comparison_grp=3-dose;
	%let suffix=4dose_3dose;
	%let vars_num=time_interval_dose3;
	* BA.2;
	%let vars_class=/*hist_covid_diag*/ hist_covid_test VA_grp_2grps /*medical_center_and_area*/ index_month; *Collapsed VA_grp further;
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA2, variant_var=variant_subgroup_sgtf, condition_model=0); 
	* BA.4/5;
	%let vars_class=/*hist_covid_diag*/ hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month; 
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA4BA5, variant_var=variant_subgroup_sgtf, condition_model=0);
	* Combine results;
    %TND_combine(table9_&suffix);

/*** 4-dose vs. unvac: VE against covid_hosp ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=4-dose;
	%let comparison_grp=Unvaccinated;
	%let suffix=4dose_unvac;
	%let vars_num=;
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month_MarApr; *Collapsed index_month; 
	%TND(hospitalsets,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA2, variant_var=variant_subgroup_sgtf, condition_model=0); 
	* BA.4/5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month; 
	%TND(hospitalsets2,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA4BA5, variant_var=variant_subgroup_sgtf, condition_model=0); *Collapsed race_eth in hospitalsets2;
	* Combine results;
    %TND_combine(table9_&suffix);

**********************************
** Table 10-11. Sensitivity analysis without IC;

data subgroupflag_nonIC;
	set subgroupflag;
	where IC=0;
run;
data hospitalsets_nonIC;
	set hospitalsets;
	where IC=0;
	if index_month in ('2022-02','2022-03','2022-04','2022-05') then index_month_FebtoMay='Feb-May'; else index_month_FebtoMay=index_month;
run;
data hospitalsets_nonIC2; set hospitalsets_nonIC;
	if race_eth in ('Non-Hispanic Asian','Non-Hispanic Black','Other/Unknown') then race_eth='Other/Unknown';
run;
data hospitalsets_nonIC3; set hospitalsets_nonIC2;
	if age_grp in ('18-44 yo','45-64 yo') then age_grp='18-64 yo';
run;
proc freq data=hospitalsets_nonIC3;
	table age_grp race_eth;
run;

** Non-IC sensitivity analyses for covid_dx;

/*** 3-dose vs. unvac: VE against covid_dx ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=3-dose;
	%let comparison_grp=Unvaccinated;
	%let suffix=3dose_unvac;
	%let vars_num=;
	* BA.1; 
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month; 
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA1, 	condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA1, 	condition_model=0, time_since_vac=1); 
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2, 	condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2,	condition_model=0, time_since_vac=1);
	* BA.2.12.1;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2121, condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2121, condition_model=0, time_since_vac=1);
	* BA.4;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month; 
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4, 	condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4, 	condition_model=0, time_since_vac=1);
	* BA.5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA5, 	condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA5, 	condition_model=0, time_since_vac=1);
	* Combine results;
    %TND_combine(table10_&suffix);

/*** 3-dose vs. 2-dose: VE against covid_dx  ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=3-dose;
	%let comparison_grp=2-dose;
	%let suffix=3dose_2dose;
	%let vars_num=time_interval_dose2;
	* BA.1; 
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA1, 	condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA1, 	condition_model=0, time_since_vac=1); 
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2, 	condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2,	condition_model=0, time_since_vac=1);
	* BA.2.12.1;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2121, condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2121, condition_model=0, time_since_vac=1);
	* BA.4;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month; 
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4, 	condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4, 	condition_model=0, time_since_vac=1);
	* BA.5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month;  
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA5, 	condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA5, 	condition_model=0, time_since_vac=1);
	* Combine results;
    %TND_combine(table10_&suffix);

/*** 4-dose vs. 3-dose: VE against covid_dx  ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=4-dose;
	%let comparison_grp=3-dose;
	%let suffix=4dose_3dose;
	%let vars_num=time_interval_dose3;
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2, 	condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2,	condition_model=0, time_since_vac=1);
	* BA.2.12.1;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2121, condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2121, condition_model=0, time_since_vac=1);
	* BA.4;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month; 
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4, 	condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4, 	condition_model=0, time_since_vac=1);
	* BA.5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month; 
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA5, 	condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA5, 	condition_model=0, time_since_vac=1);
	* Combine results;
    %TND_combine(table10_&suffix);

/*** 4-dose vs. unvac: VE against covid_dx  ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=4-dose;
	%let comparison_grp=Unvaccinated;
	%let suffix=4dose_unvac;
	%let vars_num=;
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp medical_center_and_area index_month; 
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2, 	condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2,	condition_model=0, time_since_vac=1);
	* BA.2.12.1;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month_MarApr; *Collapsed index_month;
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2121, condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA2121, condition_model=0, time_since_vac=1);
	* BA.4;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month_AprMay; *Collapsed index_month;
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4, 	condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA4, 	condition_model=0, time_since_vac=1);
	* BA.5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp /*medical_center_and_area*/ index_month; 
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA5, 	condition_model=0);
	%TND(subgroupflag_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_dx, BA5, 	condition_model=0, time_since_vac=1);
	* Combine results;
    %TND_combine(table10_&suffix);


** Non-IC sensitivity analyses for covid_hosp;

/*** 3-dose vs. unvac: VE against covid_hosp ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=3-dose;
	%let comparison_grp=Unvaccinated;
	%let suffix=3dose_unvac;
	%let vars_num=;
	* BA.1; 
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month_MarAprMay; *Collapsed index_month;
	%TND(hospitalsets_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA1, condition_model=0); 
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month; 
	%TND(hospitalsets_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA2, condition_model=0);
	* BA.4/5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month; 
	%TND(hospitalsets_nonIC2,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA4BA5, condition_model=0); *Collapsed race/eth;
	* Combine results;
    %TND_combine(table11_&suffix);

/*** 3-dose vs. 2-dose: VE against covid_hosp ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=3-dose;
	%let comparison_grp=2-dose;
	%let suffix=3dose_2dose;
	%let vars_num=time_interval_dose2;
	* BA.1; 
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month_FebtoMay; *Collapsed index_month;
	%TND(hospitalsets_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA1, condition_model=0); 
	* BA.2;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month; 
	%TND(hospitalsets_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA2, condition_model=0);
	* BA.4/5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month; 
	%TND(hospitalsets_nonIC2,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA4BA5, condition_model=0); *Collapsed race/eth;
	* Combine results;
    %TND_combine(table11_&suffix);

/*** 4-dose vs. 3-dose: VE against covid_hosp ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=4-dose;
	%let comparison_grp=3-dose;
	%let suffix=4dose_3dose;
	%let vars_num=time_interval_dose3;
	* BA.2;
	%let vars_class=/*hist_covid_diag*/ hist_covid_test VA_grp_2grps /*medical_center_and_area*/ index_month; *Collapsed VA_grp further;
	%TND(hospitalsets_nonIC,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA2, condition_model=0); 
	* BA.4/5;
	%let vars_class=/*hist_covid_diag*/ hist_covid_test VA_grp_3grps /*medical_center_and_area*/ /*index_month*/; *Dropped index month (only 2 grps, still can't converge); 
	%TND(hospitalsets_nonIC3,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA4BA5, condition_model=0); *Collapsed race/eth and age_grp;
	* Combine results;
    %TND_combine(table11_&suffix);

/*** 4-dose vs. unvac: VE against covid_hosp ***/
	* Clear any results tables;
    proc datasets lib=work; delete results:; run; 
	* Set macro values;
	%let exposed_grp=4-dose;
	%let comparison_grp=Unvaccinated;
	%let suffix=4dose_unvac;
	%let vars_num=;
	* BA.2;
	%let vars_class=/*hist_covid_diag*/ hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month_MarApr; *Collapsed index_month;
	%TND(hospitalsets_nonIC2,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA2, condition_model=0); 
	* BA.4/5;
	%let vars_class=hist_covid_diag hist_covid_test VA_grp_3grps /*medical_center_and_area*/ index_month; 
	%TND(hospitalsets_nonIC2,&exposed_grp,&comparison_grp,&suffix,covid_hosp, BA4BA5, condition_model=0); *Collapsed race_eth in hospitalsets2;
	* Combine results;
    %TND_combine(table11_&suffix);
