function Solution = LocalSearch(Pop,D,M)
N = length(Pop);
Solution = [];
for i=1:N
    Parent = Pop(i).domains;
    Solution = [Solution Parent];
end
Solution = GA(Solution);
end