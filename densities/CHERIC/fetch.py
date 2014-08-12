import urllib

for i in range(2, 2000):
    path = "http://www.cheric.org/research/kdb/hcprop/showprop.php?cmpid=%d" % i
    filename = "./pages/page%d.html" % i
    urllib.urlretrieve(path, filename)
