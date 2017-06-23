"""Routines to make input and outpu for unit tests."""

"""Numerical methods for mathematical functions neede for the program"""

#This method finds the factorial of a given number(num).
#The method calls for an integer input num
#The method returns the factorial(represented as fact in the method)
def factorial(num):
    """Finds the factorial of the input integer.
    
    :arg num: an integer
    """
    #If the number provided is zero then the factorial is 1
    if num == 0:
	fact = 1
    #Otherwise set fact to 1 and begin finding the factorial r is
    #used to find each num-n for n=0 to n=num each value in r is
    #then used to compute the factorial
    else:
	fact = 1
	r = range(1,num+1)
	for i in r:
	    fact *= i
    return fact

#This method finds a binomial coefficient that is it calculates n
#choose r The method calls for two integer values n and r and the
#method factorial The method returns the binomial factorial as result
def binomial_coefficient(n,r):
    """Finds the binomial coefficient, n choose r, for a given set of
    integers.
    
    :arg n: an integer
    :arg r: an integer
    """
    #If r is less than zero then the binomial coefficient is zero
    if r < 0:
	result = 0
    #Otherwise use the factorial method to calculate the binomial
    #coefficient n!/r!/(n-r)!
    else:
	result = factorial(n) / factorial(r) / factorial(n-r)
    return result


def group(gen):
    """Generates an entire group using the specified generators by applying generators
    to each other recursively.

    :arg gen: a list of generators as integers.
    """
    from operator import itemgetter as iget
    def g_apply(operations, source, groupi=None):
        """Applies the specified group operations to the source list of elements and then
        appends it to the group if it is unique.
        """
        result = list(iget(*operations)(source))
        if groupi is not None and result not in groupi:
            groupi.append(result)
        return result

    #Make sure the group is zero-based for python.
    if not 0 in gen[0]:
        ngens = [list(map(lambda e: e-1, g)) for g in gen]
    else:
        ngens = gen

    groupi = []
    for i in ngens:
        for j in ngens: #filter(lambda k: k!=i, ngens):
            c = g_apply(i, j, groupi)
            d = g_apply(i, c, groupi)
            while d != c:
                d = g_apply(i, d, groupi)

    while True:
        group2 = []
        for i in ngens:
            for h in groupi:
                d = g_apply(i, h)
                if d not in groupi and d not in group2:
                    group2.append(d)

        groupi.extend(group2)
        if len(group2) == 0:
            break
    return(groupi)

cases = ["120a", "120b", "Z10", "1110", "576_9", "228_8", "216_12",
         "1152_12", "256_16", "16_8", "64_8", "360_6", "20_10", "50_10", "120_10", "32_16",
         "1024_16", "441_20", "40_20", "200_20", "2500_20", "118_44_fcc", "118_44_32_fcc",
         "111_22_fcc", "111_13_fcc", "111_53_fcc", "111_233_fcc", "111_2222_fcc", "111_512_fcc",
         "651_88_fcc", "651_448_fcc", "651_646_fcc", "651_4444_fcc", "114_22_sc", "118_26_sc",
         "118_44_sc", "1110_55_sc", "1110_433_sc", "1116_88_sc", "1116_655_sc", "1116_4426_sc",
         "114_22_bcc", "118_44_bcc", "118_323_bcc", "1116_88_bcc", "113_33_hcp", "113_222_hcp",
         "116_66_hcp", "116_633_hcp", "116_3333_hcp"]

from os import system, chdir

for case in cases:
    print case

    system("cp -r ~/Documents/research/new_uncle/enum4/unittests/tree_class.initializeTree.g49/tests/initializeTree.{0}/initializeTree.out ~/Documents/research/new_uncle/enum4/tests/tree_class/initializeTree.out.{0}".format(case))
    # with open(case,"r") as cf:
    #     title = cf.readline()
    #     k = cf.readline()
    #     colors = [int(i) for i in cf.readline().split()]
    #     ng = cf.readline()
    #     gens = []
    #     for i in range(int(ng)):
    #         gens.append([int(j) for j in cf.readline().split()])

    # generators = [[j-1 for j in i] for i in gens]

    # groups = group(generators)

    # system("mkdir initializeTree.out.{}".format(case))
    # chdir("initializeTree.out.{}".format(case))
    # with open("_-k","w+") as kf:
    #     kf.write(k)
    # with open("_-n","w+") as nf:
    #     nf.write(str(sum(colors)))
    # with open("_-colors","w+") as cf:
    #     for i in colors:
    #         cf.write(str(i) + " ")
    # with open("_-loc","w+") as lf:
    #     for i in range(int(k)):
    #         lf.write("-1 ")
    # with open("_-base","w+") as bf:
    #     for i in range(int(binomial_coefficient(sum(colors),colors[0]))):
    #         bf.write("0 ")
    # with open("_-Gsize","w+") as gsf:
    #     gsf.write(str(len(groups)) +" "+ str(len(groups))+ " ")
    #     for i in range(int(k)-2):
    #         gsf.write("0 ")
    # with open("_-branches","w+") as bf:
    #     for i in range(int(k)):
    #         bf.write(str(binomial_coefficient(sum(colors[i:]),colors[i]))+ " ")
    # with open("_-G-layer-1-perms","w+") as gf:
    #     for i in groups:
    #         for j in i:
    #             gf.write(str(j+1) + " ")
    #         gf.write("\n")
            
    # chdir("../")
