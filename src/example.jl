using FormFactor

print("n1 l1 m1 n2 l2  form factor \n")
print("----------------------------------------\n")
for n1 in 1:10
    for l1 in 0:n1-1
        for m1 in 0:l1
            print(" $n1 "," $l1 "," $m1 ", " $n1 "," $l1  ",mform(n1,l1,m1,n1,l1,1.0,0.66),"\n")
        end
    end
end
