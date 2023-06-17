race = {checkEquality => false} >> o -> (fun1,fun2,args,numIterations) -> (
    time M = (fun1 args);
    time N = (fun2 args);
    if o.checkEquality then print(M==N);
    L = for i to numIterations list (
	    M = timing fun1 args;
	    N = timing fun2 args;
	    {M#0, N#0,M#0>N#0}
	);
    Mtime = L / (i -> i#0);
    Ntime = L / (i -> i#1);
    TruesFalses = L / (i->i#2);
    T = tally TruesFalses;
    Tcat = tally (TruesFalses|{true,false});
    -- Concattenating with true false so that 
    if Tcat#false>Tcat#true then (
	print(toString(fun1) | " was faster than " | toString(fun2))
	)
    else print(toString(fun2) | " was faster than " | toString(fun1));
    Mmean= mean(Mtime);
    Nmean = mean(Ntime);
    print(toString(fun1) | " took " | toString(Mmean) | " seconds on average");
    print(toString(fun2) | " took " | toString(Nmean) | " seconds on average");
    L
    )
	    

mean = L -> (sum L)/#L
