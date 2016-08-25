import os

# cases = ["120a", "120b", "Z10", "1110", "576_9", "228_8", "216_12", "1152_12", "256_16", "16_8", "64_8", "360_6", "20_10", "50_10", "120_10", "32_16", "1024_16", "441_20", "40_20", "200_20", "2500_20", "118_44_fcc", "118_44_32_fcc", "111_22_fcc", "111_13_fcc", "111_53_fcc", "111_233_fcc", "111_2222_fcc", "111_512_fcc", "651_88_fcc", "651_448_fcc", "651_646_fcc", "651_4444_fcc", "114_22_sc", "118_26_sc", "118_44_sc", "1110_55_sc", "1110_433_sc", "1116_88_sc", "1116_655_sc", "1116_4426_sc", "114_22_bcc", "116_633_hcp", "116_3333_hcp", "118_44_bcc", "118_323_bcc", "1116_88_bcc", "113_33_hcp", "113_222_hcp", "116_66_hcp"]

cases = range(1,23)

for case in cases:
    print("c",case)
    layer = 1
    growing = True
    os.chdir("incr_loc.in.{}".format(case))
    while growing:
        if os.path.isfile("_-G-layer-{}-perms".format(layer)):
            layer += 1
        else:
            growing = False

    for l in range(1,layer):
        group_size = 0
        with open("_-G-layer-{}-perms".format(l),"r") as gf:
            for line in gf:
                group_size +=1

        with open("_-A-layer-{}-perms".format(l),"w+") as af:
            for i in range(group_size):
                af.write("0\n")
    print("l",layer)
    layer2 = 0 + layer
    growing = True
    while growing:
        if os.path.isfile("_-G-layer-{}-.fpy.blank".format(layer2)):
            layer2 += 1
        else:
            growing = False
    print("l2",layer2)
    for l in range(layer,layer2):
        with open("_-A-layer-{}-.fpy.blank".format(l),"w+") as af:
            af.write('# <fortpy version="1" template="logical"></fortpy>\n F')

    with open("_-color_map".format(case),"w+") as ocf:
        ocf.write("0 0 \n 0 0")
    with open("_-narrows".format(case),"w+") as naf:
        naf.write("0")
    os.chdir("../")
        
    layer = 1
    growing = True
    os.chdir("incr_loc.out.{}".format(case))
    while growing:
        if os.path.isfile("_-G-layer-{}-perms".format(layer)):
            layer += 1
        else:
            growing = False

    for l in range(1,layer):
        group_size = 0
        with open("_-G-layer-{}-perms".format(l),"r") as gf:
            for line in gf:
                group_size +=1

        with open("_-A-layer-{}-perms".format(l),"w+") as af:
            for i in range(group_size):
                af.write("0\n")

    layer2 = 0 + layer
    growing = True
    while growing:
        if os.path.isfile("_-G-layer-{}-.fpy.blank".format(layer2)):
            layer2 += 1
        else:
            growing = False
    for l in range(layer,layer2):
        with open("_-A-layer-{}-.fpy.blank".format(l),"w+") as af:
            af.write('# <fortpy version="1" template="logical"></fortpy>\n F')

    with open("_-color_map".format(case),"w+") as ocf:
        ocf.write("0 0 \n 0 0")
    with open("_-narrows".format(case),"w+") as naf:
        naf.write("0")
    os.chdir("../")
        
