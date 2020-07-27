?- lib(ic).
?- lib(random).

?- assert(verbose).

mycount(Goal,Ct) :-
    setval(mycount,0),
    ( call(Goal), incval(mycount), fail ; 
      getval(mycount,Ct)
    ).

gen(NSF,N,File,Ct) :-
       write_gen(NSF,N,File,Ct),
       initprob(NSF,N,Ct,Cost),
       write_neighbours(N,Cost,File),
       (repeat,
        getval(limct,Curr),
        (verbose -> writeln(it(Curr)-start) ; true), 
        rand(Cost,List),
       ( NSF= nsf -> List = [X|_], X=N ; % To constrain the values so the initial domain is sufficient
                      true
        ),
        (verbose -> writeln(it(Curr)-gen-index) ; true),    
        gen_random(List,MaxVal,MinVal),
        write_iter(File,Curr,MaxVal,MinVal),  
        write_cost(Cost,File),
%        write_trans(Curr,MaxVal,Cost,File),
        incval(limct),
        getval(limct,Ct)
       ), !.

write_gen(NSF,N,File,Ct) :-
    open(File,append,S),
    nl(S),
    writeln(S,type:NSF - variable_count:N - iterations:Ct),
    close(S).


write_iter(File,Curr,MaxVal,MinVal) :-
        open(File,append,S),
        nl(S),
        writeln(S, it(Curr)-gen-random-max-min(MaxVal,MinVal)),  
        close(S).

write_neighbours(N,Cost,File) :-
    open(File,append,S),
    nl(S),
    writeln(S, new_fitness_function(N)),
    writeln(S, "For each state, a list of its neighbours"),
    dim(Cost,[Dim]),
    (for(J,1,Dim), param(S) do 
        (findall(K,n(J,K),Neighbours), writeln(S, J-Neighbours))
    ),
    close(S).

write_cost(Cost,File) :-
    open(File,append,S),
    nl(S),
    writeln(S,"For each state, its fitness value"),
    dim(Cost,[Dim]),
    (for(J,1,Dim-1), param(Cost,S) do 
        (X is Cost[J], write(S,J-X),write(S,', '))
    ),
    X is Cost[Dim],
    writeln(S,Dim-X),
    close(S).

write_trans(Curr,MaxVal,Cost,File) :-
    open(File,append,S),
    nl(S),
%    writeln(S,"Transition probabilities: For each solution, a list of neighbours, and the number of neighbours."),
    dim(Cost,[Dim]),
    write(S,it(Curr)),write(S,';solutions '),writeln(S,Dim),  
    (for(Sol,1,Dim), param(MaxVal,Cost,Dim,S) do
      C is Cost[Sol],
      ( C=:=MaxVal -> Neighbours = []
        ;
        n(Sol,K), Cost[K]>C -> 
           findall(Neigh,(n(Sol,Neigh),Cost[Neigh]>C), Neighbours) 
        ;
           findall(Neigh,(n(Sol,Neigh),Cost[Neigh]=:=C), Neighbours) 
      ),
      find_fract(Dim,Neighbours,NCList),
      write_aldeida(S,NCList)
%      find_error(Cost,NCList)
%        write(S, J-NCList),write(S,', ')
    ),
    nl(S),
    close(S).

find_error(Cost,NCList) :-
   nth(N,NCList,Neighbours),
   Cost[N] = C,
   nth(M,Neighbours,Fract), Fract>=0.1,
   Cost[M] < C,
   writeln(error(N,M)), !,
   abort;
   writeln('NCList-OK').

nth(1,[H|_],H).
nth(N,[_|T],V) :-
   nth(M,T,V), N is M+1.

write_aldeida(S,NCList) :-
    ( foreach(_K-Fract,NCList),param(S) do write(S,Fract),write(S,' ') ),
    nl(S).

find_fract(Dim,Neighbours,NCList) :-
        length(Neighbours,Len),
        (Len = 0 -> Fract = 0.0 ; 
                    IntFrac is 1000 // Len, Fract is IntFrac/1000
        ),
        (for(K,1,Dim), foreach(K-Fraction,NCList), param(Neighbours,Fract)
        do
           (memberchk(K,Neighbours) -> Fraction is Fract ; Fraction is 0.0)
        ).

rand(Cost,RList) :-
    dim(Cost,[Dim]),
    Cost =.. [[]|List],
    recr(Dim,List,RList).

recr(0,[],[]) :- !.
recr(N,List,[H|STail]) :-
    random(N,X),
    apx(X,List,Rest,H),
    M is N-1,
    recr(M,Rest,STail).

apx(0,[H|T],T,H) :- !.
apx(X,[H|T],[H|R],V) :-
    Y is X-1,
    apx(Y,T,R,V).

initprob(NSF,N,Ct,Cost) :-
   (var(Ct) -> writeln(no_count(Ct)), abort ; true),
   Dim is 2^N,
   dim(Cost,[Dim]),
   ( NSF= nsf -> Cost :: 0..2*N;  % Enough to allow maximum difference from any initial point to be N
                 Cost::0..N
   ),
   gen_nlt(N),
   (NSF=nsf -> nc(Cost) ; NSF=no_nsf -> true ; writeln(wrong(NSF)), abort ),
   setval(limct,0).

gen_random(List,MaxVal,MinVal) :-
    (foreach(X,List),fromto(0,This,Next,MaxVal),fromto(1000,TMin,NMin,MinVal)
    do
       indomain(X,random),
       (X>This -> Next=X ; Next=This),
       (X<TMin -> NMin=X ; NMin=TMin)
    ), !.
%    search(List,0,input_order,indomain_random,complete,[]), !.

% Neighbour Constraiint

nc(Cost) :-
   findall([X,Y],n(X,Y),Neighbours),
   (foreach([N1,N2],Neighbours), param(Cost)  
        do Cost[N2] - Cost[N1] #=< 1, Cost[N1]-Cost[N2] #=< 1
   ).


gen_nlt(N) :-
   assert(nlt(0,0)),
   retractall(nlt(_,_)),
   (for(I,1,N), param(N) do gen_nlt(I,N)).

n(X,Y) :- nlt(X,Y) ; nlt(Y,X).
gen_nlt(I,N) :-
   dim(X1,[N]),dim(X2,[N]), 
   X1 :: 0..1, X2 :: 0..1,
   X1[I] #= 0, X2[I] #= 1,
   (for(J,1,N), param(I,X1,X2) do 
        (J=I -> true ; X1[J] #= X2[J])
   ),
   ( labeling(X1),
     mycalc(N,X1,V1), mycalc(N,X2,V2),
     assert(nlt(V1,V2)),
     fail ;
     true
   ).

mycalc(N,Array,Val) :-
   for(I,1,N), fromto(1,TV,NV,Val), param(Array) do
        NV is TV + Array[I]*(2^(I-1)).




