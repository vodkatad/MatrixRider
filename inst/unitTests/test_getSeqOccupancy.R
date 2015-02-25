# library(RUnit)
# source('inst/unitTests/test_getSeqOccupancy.R')


# I do not like that much to depend on JASPAR2014, maybe in the future build ad
# hoc matrixes to have more robust/reliable tests.
test_getSeqOccupancy_singlePWM_manyCutoffs <- function() {
    library(JASPAR2014)
    library(TFBSTools)
    library(Biostrings)
    pfm <- getMatrixByID(JASPAR2014,"MA0004.1")
    sequence <- DNAString("CACGTG")
    cutoffs <- seq(0,1,0.1)
    given <- vapply(cutoffs, FUN=function(x) {
        getSeqOccupancy(sequence=sequence, pfm=pfm, cutoff=x)
    }, FUN.VALUE=1)
    wanted <- rep(1470.946,11)
    checkEqualsNumeric(given, wanted, tolerance=1e-5)
}

# More "standalone" tests that does not depend on JASPAR2014.
test_getSeqOccupancy_singlePWM_manyCutoffs_manualMat <- function() {
    library(TFBSTools)
    library(Biostrings)
    pfm <- PFMatrix(ID="Unknown", name="Unknown", matrixClass="Unknown",
                    strand="*", bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                    tags=list(), profileMatrix=matrix(c(4L,  19L, 0L,  0L,  0L,
                            0L,
                            16L, 0L,  20L, 0L,  0L,  0L,
                            0L,  1L,  0L,  20L, 0L,  20L,
                            0L,  0L,  0L,  0L,  20L, 0L),
                            byrow=TRUE, nrow=4,
                            dimnames=list(c("A", "C", "G", "T"))))
    sequence <- DNAString("CACGTG")
    cutoffs <- seq(0,1,0.1)
    given <- vapply(cutoffs, FUN=function(x) {
        getSeqOccupancy(sequence=sequence, pfm=pfm, cutoff=x)
    }, FUN.VALUE=1)
    wanted <- rep(1470.946,11)
    checkEqualsNumeric(given, wanted, tolerance=1e-5)
}

# is there a way to load the needed packages only once? Maybe they could be but 
# in runTests.R?
test_getSeqOccupancy_manyPWM_totAff <- function() {
    library(JASPAR2014)
    library(TFBSTools)
    library(Biostrings)
    opts<-list()
    opts['collection']='CORE'
    opts['matrixtype']='PFM'
    opts['all_versions']=FALSE
    opts['tax_group']='vertebrates'
    mat <- getMatrixSet(JASPAR2014, opts)
    sequence <- DNAString("CACGTGCTGAAATGTCATGCATGCTAGCGTTTTGCCATGTCGATCTGACTCGTAGTGCTGTCGTTTGCGGATCGTTCGAGCTGCTGCGTATCGTGC")
    given <- getSeqOccupancy(sequence=sequence, pfm=mat, cutoff=0)
    wanted <- structure(c(1507.81772494717, 216.555152004053, 5.1231094462624,
                        143.266000843078, 44.3141298556815, 
                        23.1218205263211, 119.19900616835,
                        138.885871138134, 5.30424036430996,
                        24.6533423461804, 28.1719256935498,
                        134.151580040822, 61.9866353334519,
                        26.920024332854, 6.67084700228316,
                        5.33615336765674, 3.94579032756123,
                        133.822499797666, 1.14085340664733,
                        1869.13667373677, 710.524842986656,
                        51.9198177838886, 26.0874196264541,
                        23.6119262113813, 24.4418595163519,
                        0.383750678036033, 786.356042971725,
                        8.65920513177419, 41.407546352327, 
                        11.1899919743433, 29.7952444309024,
                        0.243418724077108, 115.454674552202,
                        7.7764295829428, 2.91011542069716,
                        15.6189770467564, 38.3339243562868,
                        66.2228632111535, 26.6224427573005,
                        67.8233502745905, 60.8114917203315,
                        736.786150813535, 128.813111240648,
                        5.05292357982437, 112.047772066333,
                        31.2516651070403, 12.4685697428492,
                        0.725199686006218, 93.2181781794895,
                        56.2519636851859, 7.86572499271935,
                        0.832064754971665, 390.630066392318,
                        277.764626599316, 50.091577717651,
                        22.1533393823669, 2.70023528612808,
                        116.452364957258, 60.4564942066966,
                        21.1588217136918, 316.352067469659,
                        0.672802854323599, 32.9498817882305,
                        0.0522768271356331, 51.9432783559069,
                        7.1507226053144e-07, 4.86368233329524,
                        0.0465332631131116, 0.000910248596604364,
                        17.4458294041035, 0.103539007105645,
                        1253.49209002998, 7.78271122775548,
                        14.6338991181658, 210.977975917038,
                        21.6381482021818, 5.47644468994859,
                        97.783199023199, 48.7627819601915,
                        273.385054968354, 1.75170145946744,
                        89.7720965309201, 99.4724121159676,
                        0.241595909432787, 53.6159955989287,
                        57.6467962321399, 195.745863362006,
                        21.1169351961258, 61.0508095447226,
                        13.6338526893635, 4.40703787815226,
                        0.0605322677528101, 0.585981182633185,
                        2.24450118002449, 2.98652178879596,
                        0.000175847458049154, 0.00164688340497543,
                        0.104175907337801,
                        0.00886930757238849, 0.000340453185600774,
                        0.00691921268835147,
                        0.177049538732553, 0.0347323378515403, 
                        0.00413472726806136, 0.0415440596883874,
                        0.00550602669636055, 0.135104348594409,
                        0.832649991913259, 1.2012942363265,
                        0.871565811433871, 0.00227214329349484, 
                        0.239544812642387, 0.346685428899121,
                        0.0154132708398723, 24.2786906350577,
                        6.37097890769897, 0.0458158964051514,
                        90.8239163831229, 0.017258972229131, 
                        0.545536114597832, 0.487440858989092,
                        0.284521503649544, 0.000281958921606087,
                        0.015886072707717, 0.919121714754885,
                        17.9475193227739, 16.9659484155644,
                        0.153458694453079, 32.3100043441455,
                        10.4007346402548, 1.31227227347278,
                        7.01487767100741, 0.000719727015688071,
                        0.0559454565831749, 0.136608906795682,
                        0.0059914735105239, 0.0249386077339904,
                        0.174508509236544, 0.0139069904484719,
                        5.7258580721095, 0.0542235603706622,
                        66.4745215979321, 0.414822527016949,
                        6.02987553489201, 0.0286658510344524,
                        0.336298941636761, 0.000944701685256303,
                        0.129743605927054, 0.0100255277124301,
                        6.7853006317582e-05, 0.000130327567499375,
                        1.744739458302, 1.0321531780572,
                        0.00528610768330774, 0.00212644610951357, 
                        2.9804404754525, 0.000182906762570604,
                        0.0908920526837976, 0.00267757081520737,
                        16.6801820154201, 0.456976467490772,
                        0.979301759925958, 2.871329053895e-06,
                        0.000635725514267653,
                        0.669611458102948, 477.286815851959,
                        0.0281962977251376, 0.000599561278291718,
                        0.00363795475349825, 0.0211697177941742,
                        0.000469517815982694,
                        0.458235794752639, 0.0204842332950218, 
                        0.186336846766483, 0.00456749366423416,
                        0.126414804194817, 2.28543644615045,
                        4091.890809206, 12.5607977988511,
                        0.00874606664842197, 0.0948719721322671,
                        316.295587994671, 0.000446809852082557,
                        0.00331331809611546, 1.11441529486672e-06,
                        0.0252051436135396,
                        0.000211382216757059, 0.00869801591910128
                        , 0.0902657308710231,
                        0.00608590192951825, 0.103687932106769,
                        0.0165060021014958, 87.1522042887476,
                        0.00176943454809679, 11.7472273767883, 
                        11.1859220274905, 0.00491005489690446,
                        0.455107293092024, 1.47660412710369,
                        3.4199837849597, 103.906635699472,
                        0.00793317609363519, 4.53205749009618e-07,
                        5.13518069769366,
                        0.209700495324655), 
                        .Names = c("MA0004.1", "MA0006.1", "MA0009.1",
                                   "MA0017.1", "MA0019.1", "MA0025.1", 
                                   "MA0027.1", "MA0028.1", "MA0029.1",
                                   "MA0030.1", "MA0031.1", "MA0032.1", 
                                   "MA0033.1", "MA0038.1", "MA0040.1",
                                   "MA0041.1", "MA0042.1", "MA0043.1",
                                   "MA0046.1", "MA0048.1", "MA0051.1",
                                   "MA0056.1", "MA0057.1", "MA0059.1",
                                   "MA0063.1", "MA0066.1", "MA0067.1",
                                   "MA0068.1", "MA0069.1", "MA0070.1",
                                   "MA0071.1", "MA0072.1", "MA0073.1",
                                   "MA0074.1", "MA0075.1", "MA0077.1",
                                   "MA0078.1", "MA0081.1", "MA0084.1",
                                   "MA0087.1", "MA0088.1", "MA0089.1",
                                   "MA0090.1", "MA0091.1", "MA0092.1",
                                   "MA0101.1", "MA0107.1", "MA0108.2",
                                   "MA0109.1", "MA0111.1", "MA0115.1",
                                   "MA0116.1", "MA0117.1", "MA0119.1", 
                                   "MA0122.1", "MA0124.1", "MA0125.1",
                                   "MA0130.1", "MA0131.1", "MA0132.1",
                                   "MA0133.1", "MA0135.1", "MA0136.1",
                                   "MA0139.1", "MA0142.1", "MA0149.1",
                                   "MA0062.2", "MA0039.2", "MA0138.2",
                                   "MA0002.2", "MA0047.2", "MA0112.2",
                                   "MA0065.2", "MA0151.1", "MA0152.1",
                                   "MA0153.1", "MA0155.1", "MA0156.1",
                                   "MA0157.1", "MA0158.1", "MA0159.1",
                                   "MA0160.1", "MA0161.1", "MA0163.1",
                                   "MA0164.1", "MA0018.2", "MA0099.2",
                                   "MA0259.1", "MA0442.1", "MA0141.2",
                                   "MA0145.2", "MA0146.2", "MA0461.1",
                                   "MA0462.1", "MA0463.1", "MA0464.1", 
                                   "MA0465.1", "MA0466.1", "MA0467.1",
                                   "MA0468.1", "MA0469.1", "MA0470.1",
                                   "MA0471.1", "MA0472.1", "MA0473.1",
                                   "MA0474.1", "MA0475.1", "MA0476.1",
                                   "MA0477.1", "MA0478.1", "MA0479.1",
                                   "MA0480.1", "MA0481.1", "MA0482.1", 
                                   "MA0483.1", "MA0484.1", "MA0485.1",
                                   "MA0486.1", "MA0488.1", "MA0489.1", 
                                   "MA0490.1", "MA0491.1", "MA0492.1",
                                   "MA0493.1", "MA0494.1", "MA0495.1", 
                                   "MA0496.1", "MA0497.1", "MA0498.1",
                                   "MA0499.1", "MA0500.1", "MA0501.1", 
                                   "MA0502.1", "MA0503.1", "MA0504.1",
                                   "MA0505.1", "MA0506.1", "MA0507.1", 
                                   "MA0508.1", "MA0509.1", "MA0510.1",
                                   "MA0511.1", "MA0512.1", "MA0513.1",
                                   "MA0514.1", "MA0515.1", "MA0516.1",
                                   "MA0517.1", "MA0518.1", "MA0519.1", 
                                   "MA0520.1", "MA0521.1", "MA0522.1",
                                   "MA0523.1", "MA0524.1", "MA0525.1", 
                                   "MA0526.1", "MA0527.1", "MA0528.1",
                                   "MA0007.2", "MA0102.3", "MA0024.2", 
                                   "MA0154.2", "MA0162.2", "MA0076.2",
                                   "MA0258.2", "MA0098.2", "MA0148.3", 
                                   "MA0035.3", "MA0036.2", "MA0037.2",
                                   "MA0114.2", "MA0050.2", "MA0058.2",
                                   "MA0052.2", "MA0100.2", "MA0147.2",
                                   "MA0104.3", "MA0150.2", "MA0105.3", 
                                   "MA0060.2", "MA0014.2", "MA0080.3",
                                   "MA0143.3", "MA0079.3", "MA0083.2", 
                                   "MA0137.3", "MA0144.2", "MA0140.2",
                                   "MA0003.2", "MA0106.2", "MA0093.2", 
                                   "MA0095.2", "MA0103.2", "MA0591.1",
                                   "MA0595.1", "MA0596.1", "MA0597.1",
                                   "MA0598.1", "MA0599.1", "MA0600.1",
                                   "MA0113.2"))
    
    checkEqualsNumeric(given, wanted, tolerance=1e-5)
}

# Trivial error-exception checking.
test_getSeqOccupancy_badCutoff <- function() {
    library(TFBSTools)
    library(Biostrings)
    pfm <- PFMatrix(ID="Unknown", name="Unknown", matrixClass="Unknown",
                    strand="*", bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                    tags=list(), profileMatrix=matrix(c(4L,  19L, 0L,  0L,  0L,  
                            0L, 16L, 0L,  20L, 0L,  0L,  0L,
                            0L,  1L,  0L,  20L, 0L,  20L,
                            0L,  0L,  0L,  0L,  20L, 0L),
                            byrow=TRUE, nrow=4,
                            dimnames=list(c("A", "C", "G", "T"))))
    sequence <- DNAString("CACGTG")
    obs <- tryCatch(getSeqOccupancy(sequence=sequence, pfm, 10), 
                    error=function(e) e)
    checkEquals(obs$message, paste0("Wrong argument to getSeqOccupancy, ",
                "'cutoff' has to be between 0 and 1! (0 <= cutoff <= 1)"))
}

# This is a test on C code but as long as it's strictly dependent on the user's 
# usage of getSeqOccupancy it could be put
# here even if it's not an UnitTest strictly speaking but an interaction one.
test_getSeqOccupancy_badDNAString <- function() {
    library(TFBSTools)
    library(Biostrings)
    pfm <- PFMatrix(ID="Unknown", name="Unknown", matrixClass="Unknown",
                    strand="*", bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                    tags=list(), profileMatrix=matrix(c(4L,  19L, 0L,  0L,  0L, 
                            0L, 16L, 0L,  20L, 0L,  0L,  0L,
                            0L,  1L,  0L,  20L, 0L,  20L,
                            0L,  0L,  0L,  0L,  20L, 0L),
                            byrow=TRUE, nrow=4,
                            dimnames=list(c("A", "C", "G", "T"))))
    sequence <- DNAString("CACSTG")
    obs <- tryCatch(getSeqOccupancy(sequence=sequence, pfm, 1), 
                    error=function(e) e)
    checkEquals(obs$message, paste0("Wrong argument to getSeqOccupancy, ", 
                "'sequence' must be ",
                "based on a    restricted alphabet with only 'A','C','G','T' ",
                "and 'N'"))
}

test_Ccode_altogether <- function() {
    res <- .Call("run_tests")
    checkEqualsNumeric(res, 0);
}