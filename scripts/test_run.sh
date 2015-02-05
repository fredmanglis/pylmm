#! /bin/bash
#  ./scripts/test_run.sh &> test.out

/usr/bin/time -v env PYTHONPATH=. python scripts/pylmmKinship.py -v --bfile data/test_snps.132k.clean.noX out.kin --test --test -t 1
/usr/bin/time -v env PYTHONPATH=. python scripts/pylmmGWAS.py -v --bfile data/snps.132k.clean.noX --kfile data/snps.132k.clean.noX.pylmm.kin --phenofile data/snps.132k.clean.noX.fake.phenos out.foo --test -t 1
md5sum *

/usr/bin/time -v env PYTHONPATH=. python scripts/pylmmKinship.py -v --bfile data/test_snps.132k.clean.noX out.kin --test --test
/usr/bin/time -v env PYTHONPATH=. python scripts/pylmmGWAS.py -v --bfile data/snps.132k.clean.noX --kfile data/snps.132k.clean.noX.pylmm.kin --phenofile data/snps.132k.clean.noX.fake.phenos out.foo --test
md5sum *

/usr/bin/time -v env PYTHONPATH=. python scripts/pylmmKinship.py -v --bfile data/test_snps.132k.clean.noX out.kin --test
/usr/bin/time -v env PYTHONPATH=. python scripts/pylmmGWAS.py -v --bfile data/snps.132k.clean.noX --kfile data/snps.132k.clean.noX.pylmm.kin --phenofile data/snps.132k.clean.noX.fake.phenos out.foo
md5sum *
