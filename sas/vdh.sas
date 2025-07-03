/* https://people.musc.edu/~bandyopd/bmtry711.11/lecture_24.pdf */

/* Heijden 2018 */
    data work.DUAL_SYSTEM;
        input id list1 $ list2 $ class1 $ class2 $ count;
        datalines;
        1 1 1 0 0 259
        2 1 1 0 1 539
        3 1 1 1 0 110
        4 1 1 1 1 177
        5 0 1 . 0 91
        6 0 1 . 1 164
        7 1 0 0 . 13898
        8 1 0 1 . 12356
        9 0 0 . . .
        ;
    run;
    
/* Heijden 2021 */
    data work.DUAL_SYSTEM;
        input id list1 $ list2 $ class1 $ class2 $ count;
        datalines;
        1  1 1 0 0 3004335
        2  1 1 0 1 31995
        3  1 1 0 . 150840
        4  1 0 0 . 38634        
        5  1 1 1 0 108189
        6  1 1 1 1 435465
        7  1 1 1 . 12405
        8  1 0 1 . 4368        
        9  1 1 . 0 16512
        10 1 1 . 1 2769
        11 1 1 . . 900
        12 1 0 . . 438        
        13 0 1 . 0 398838
        14 0 1 . 1 146976
        15 0 1 . . 24636
        16 0 0 . . .
        ;
    run;

/* tabulate data */
    title1  "Table 01: Input Data";
    proc sql;
        SELECT 
            * 
        FROM
            work.DUAL_SYSTEM;
    quit;

    title1  "Table 02: Input Data Crosstab";
    proc tabulate data = work.DUAL_SYSTEM missing;
        class
            list1 list2 / descending;
        class
            class1 class2;
        var 
            count;
        table 
            (list1*class1 ), 
            (list2*class2 )*(count=' '*sum=' '*f=COMMA16.0); 
    run;

    title1 "Table 03: Uniform distribution / Initial Expectations";
    proc sql;
        CREATE TABLE work.VDH_INIT AS
        SELECT
            ID,
            list1,
            list2,
            coalesce(ds.class1, c1.class1) as class1,
            coalesce(ds.class2, c2.class2) as class2,
            ds.count,
            (
                (CASE WHEN ds.class1 IS NULL THEN (SELECT COUNT(DISTINCT class1) FROM work.DUAL_SYSTEM) ELSE 1 END) *
                (CASE WHEN ds.class2 IS NULL THEN (SELECT COUNT(DISTINCT class2) FROM work.DUAL_SYSTEM) ELSE 1 END)
            ) AS numImputedClassesUsedBy,
            ds.count / (
                (CASE WHEN ds.class1 IS NULL THEN (SELECT COUNT(DISTINCT class1) FROM work.DUAL_SYSTEM) ELSE 1 END) *
                (CASE WHEN ds.class2 IS NULL THEN (SELECT COUNT(DISTINCT class2) FROM work.DUAL_SYSTEM) ELSE 1 END)
            ) AS count_init
        FROM
            work.DUAL_SYSTEM AS ds
        LEFT JOIN (SELECT DISTINCT class1 FROM work.DUAL_SYSTEM WHERE class1 IS NOT NULL) AS c1 ON ds.class1 IS NULL
        LEFT JOIN (SELECT DISTINCT class2 FROM work.DUAL_SYSTEM WHERE class2 IS NOT NULL) AS c2 ON ds.class2 IS NULL
        ;
    quit;

    title1  "Table 04: Uniform distribution / Initial Expectations Crosstab";
    proc tabulate data = work.VDH_INIT missing;
        class
            list1 list2 / descending;
        class
            class1 class2;
        format
            list1 list2 $LIST.;
        var 
            count_init;
        table (list1*class1 ), (list2*class2 )*(count_init=' '*sum=' '*f=COMMA16.0); */ printmiss;
    run;
    
/* genmod macro */
    %macro em(max_iter, tol);   
        %global iter;
        %global mse;
        %let iter = 1;
		%let mse  = 1;
		
        ods select none; proc genmod data = work.VDH_INIT;
            class list1 class1 list2 class2;
            model count_init = list1|class2 class1|class2 list2|class1 / dist = poisson link = log maxiter = 100;
            output out = work.VDH predicted = expectation&iter.;
        run; ods select all; 

        proc sql nowarnrecurs;
            create table work.VDH AS
            select 
                vdh.*,
                CASE
                    WHEN vdh.list1='0' AND vdh.list2='0' THEN .
                    ELSE ROUND(COALESCE( (vdh.expectation&iter. / mrg.expTot) * orgTot, vdh.expectation&iter.)) 
                    END AS expectation&iter._tmp
            from 
                work.VDH AS vdh
            LEFT JOIN (
                    SELECT ID, SUM(expectation&iter.) AS expTot, SUM(count_init) AS orgTot FROM work.VDH GROUP BY ID
                ) AS mrg 
            ON
                vdh.ID = mrg.ID                
            ;
        quit;
        
        %do %while( &iter < &max_iter and &mse > &tol);
            %let iter_prev = &iter;
            %let iter = %eval(&iter+1);
            
            ods select none; proc genmod data = work.VDH;
                class list1 class1 list2 class2;
                model expectation&iter_prev._tmp = list1|class2 class1|class2 list2|class1 / dist = poisson link = log maxiter = 100;
                output out = work.VDH predicted = expectation&iter.;
            run; ods select all;

            proc sql nowarnrecurs;
                create table work.VDH AS
                select 
                    vdh.*,
                    CASE
                        WHEN vdh.list1='0' AND vdh.list2='0' THEN .
                        ELSE ROUND(COALESCE( (vdh.expectation&iter. / mrg.expTot) * orgTot, vdh.expectation&iter.)) 
                        END AS expectation&iter._tmp,
                    (vdh.expectation&iter. / mrg.expTot) AS expectation&iter._ratio
                from 
                    work.VDH AS vdh
                LEFT JOIN (
                        SELECT ID, SUM(expectation&iter.) AS expTot, SUM(count_init) AS orgTot FROM work.VDH GROUP BY ID
                    ) AS mrg 
                ON
                    vdh.ID = mrg.ID                
                ;
            quit;
            
            proc sql noprint;
                select mean((expectation&iter. - expectation%eval(&iter.-1))**2) INTO: mse FROM work.vdh;
            quit;
        %end;
        %if &iter = &max_iter %then %do;
            %put WARNING: No convergence - Final MSE: &mse.;
        %end;
    %mend;
    
/* Run EM Algorithm */
    %em(1000, 0.0001);
    


    title1  "Table 05: Output Data Crosstab";
    proc tabulate data = work.VDH missing;
        class
            list1 list2 / descending;
        class
            class1 class2;
        format
            list1 list2 $LIST.;
        var 
            expectation&iter._tmp;
        table (list1*class1), (list2*class2)*(expectation&iter._tmp=' '*sum=' '*f=COMMA16.0); */ printmiss;
    run;

    title1  "Table 06: Output Data Final Iterations";
    proc sql;
        select
        list1, list2, class1, class2, count_init, expectation&iter., expectation&iter._tmp, expectation%eval(&iter.-1)_tmp, expectation%eval(&iter.-1), expectation&iter._ratio
        from work.vdh;
    quit;

    proc sql nowarnrecurs;
        create table work.VDH_final AS
        select 
            vdh.*,
            ROUND(COALESCE( (vdh.expectation&iter. / mrg.expTot) * orgTot, vdh.expectation&iter.)) AS final
        from 
            work.VDH AS vdh
        LEFT JOIN (
                SELECT ID, SUM(expectation&iter.) AS expTot, SUM(count_init) AS orgTot FROM work.VDH GROUP BY ID
            ) AS mrg 
        ON
            vdh.ID = mrg.ID                
        ;
    quit;

    title1  "Table 06: Final Dual System Estimate Crosstab";
        proc tabulate data = work.VDH_final missing;
            class
                list1 list2 / descending;
            class
                class1 class2;
            format
                list1 list2 $LIST.;
            var 
                final;
            table (list1*class1 all), (list2*class2 all)*(final=' '*sum=' '*f=COMMA16.0); */ printmiss;
        run;

    %put Final Iterations: &iter.;
    %put Final MSE: &mse.;