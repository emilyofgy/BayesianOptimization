/*
 Let's walk through the Latin Hypercube Sampling algorithm with a 2d example.
 
 In order to obtain a Latin hypercube sample of size T from the 2d space, 
 we devise the following algorithm: 
   1. Divide the interval of each dimension into T equally spaced subintervals. 
   2. Randomly sample two points Eta and Nu within each subinterval. 
      Call them Eta(i) and Nu(i) when they are sampled from the interval(i), 
      that is from [(i-1)/T, i/T). 
   3. At this point we have Eta(1), ..., Eta(T) and Nu(1),..., Nu(T). 
      Let i = 1 and consider Eta(i). 
      Pick an integer j at random such that 1 <= j <= T 
      and let x(i) = (Eta(i), Nu(j)) be the first point in the Latin hypercube sample. 
      Now, refrain j from further consideration, 
      increase i by 1 and repeat this step until all the Eta(i) are considered. 
      This gives us a Latin hypercube sample over the unit square that we label as x(1),... ,X(T)
*/


%macro latinHyperCubeSamp;
start latinHyperCubeSamp(seed=0, sampleSize, designTable, resultTable); 
   /* read design table */ 
   use designTable; 
   read all var {minValue maxValue} into designTable; 
   read all var {parameterName} into varNames; 
   close designTable;  
   varNames = T(varNames);
   

   /* sampling */ 
   dimension = nrow(designTable);
   Eta =j(dimension, sampleSize, .);
   do d=1 to dimension; 
      minValue = designTable[d,1];
      maxValue = designTable[d,2];
      intervals = do(minValue, maxValue, (maxValue-minValue)/sampleSize); 
      do i=1 to sampleSize;
         Eta[d, i]= intervals[1,i] + uniform(seed)*(maxValue-minValue)/sampleSize;
      end;   
   end;  
   
   lhs = j(sampleSize, dimension, .);
   mattrib lhs colname=varNames;
   selectedFlag = j(dimension, sampleSize, 0);
   do i=1 to sampleSize; 
      lhs[i, 1] = Eta[1, i];
      do d=2 to dimension;
         notSelected = loc(selectedFlag[d,] = 0);
         u = uniform(seed); 
         selection = ceil(ncol(notSelected)*u); 
         j = notSelected[1, selection];
         selectedFlag[d, j] = 1; /* refrain j from further consideration */
         Nu = Eta[d, j];
         lhs[i, d] = Nu;
      end;
   end;
   
   tbl = TableCreate(varNames, lhs); 
   call TableWriteToDataSet(tbl, resultTable);
finish; 
%mend latinHyperCubeSamp;


/*******************************************/
/* Ex.1: Sampling from 2-dimensional space */
/*******************************************/

data designTable;
   input parameterName :$32. minValue maxValue ;
   cards;
x1 0 1
x2 -2 -1
;
run;


proc iml; 
   %latinHyperCubeSamp;
   run latinHyperCubeSamp(0, 20, designTable, "lhsExampleResult1"); 
run;
quit;
 

ods graphics / width=8in height=6in;
proc sgplot data=lhsExampleResult1 cycleattrs;
   title "Latin Hypercube Example 1";
   scatter x=x1 y=x2 / markerattrs=(symbol=CircleFilled);
   xaxis grid minorgrid minor minorcount=3;
   yaxis grid minorgrid minor minorcount=3;
run; 
quit;  



/*******************************************/
/* Ex.2: Sampling from 5-dimensional space */
/*******************************************/
data designTable;
   input parameterName :$32. minValue maxValue ;
   cards;
parm1 0 1
parm2 -2 -1
parm3 10 15
parm4 0.1 0.3
parm5 100 200
;
run;


proc iml; 
   %latinHyperCubeSamp;
   run latinHyperCubeSamp(0, 10, designTable, "lhsExampleResult2"); 
run;
quit;


title "Latin Hypercube Example 2";
proc print data=lhsExampleResult2;
run;quit;  
