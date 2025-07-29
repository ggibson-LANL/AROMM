a=[1,2,3,4]
b=[2,3,6,7]

sum =[]
for i in range(len(a)):
    s = a[i]+b[i]
    sum.append(s)
    a[i]=(a[i]/s)
    b[i]=(b[i]/s)
print(sum)

print(a)
print(b)

