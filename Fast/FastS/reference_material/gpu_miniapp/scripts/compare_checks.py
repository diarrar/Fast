#!/usr/bin/env python3
"""Compare les lignes CHECK de deux sorties de la mini-app.

usage : compare_checks.py fichier1 fichier2 [--map nomA=nomB ...] [--tol 1e-12]

Pour chaque nom de noyau present dans les deux fichiers (ou mis en
correspondance par --map), compare sum, abssum et max de chaque variable
et affiche l ecart relatif. Code retour 1 si un ecart depasse la tolerance.
"""
import re
import sys


def parse(path):
    checks = {}
    rx = re.compile(r"^CHECK\s+(\S+)\s+v(\d)\s+sum=\s*(\S+)\s+abssum=\s*(\S+)\s+max=\s*(\S+)")
    for line in open(path):
        m = rx.match(line)
        if m:
            name, v = m.group(1), int(m.group(2))
            checks[(name, v)] = tuple(float(x) for x in m.group(3, 4, 5))
    return checks


def main():
    args = [a for a in sys.argv[1:]]
    tol = 1e-12
    mapping = {}
    files = []
    i = 0
    while i < len(args):
        if args[i] == "--tol":
            tol = float(args[i + 1]); i += 2
        elif args[i] == "--map":
            a, b = args[i + 1].split("="); mapping[a] = b; i += 2
        else:
            files.append(args[i]); i += 1
    if len(files) != 2:
        print(__doc__); sys.exit(2)
    c1, c2 = parse(files[0]), parse(files[1])
    worst = 0.0
    names1 = sorted({n for n, _ in c1})
    for n1 in names1:
        n2 = mapping.get(n1, n1)
        if (n2, 1) not in c2:
            continue
        for v in range(1, 6):
            a, b = c1[(n1, v)], c2[(n2, v)]
            # sum peut etre ~0 (champ conservatif) : on le compare a abssum
            rel_sum = abs(a[0] - b[0]) / max(a[1], 1e-300)
            rel_abs = abs(a[1] - b[1]) / max(a[1], 1e-300)
            rel_max = abs(a[2] - b[2]) / max(a[2], 1e-300)
            r = max(rel_sum, rel_abs, rel_max)
            worst = max(worst, r)
            flag = "" if r <= tol else "   <-- DEPASSE tol=%g" % tol
            print("%-16s -> %-16s v%d  sum %9.2e  abssum %9.2e  max %9.2e%s"
                  % (n1, n2, v, rel_sum, rel_abs, rel_max, flag))
    print("ecart relatif max : %9.2e (tol %g)" % (worst, tol))
    sys.exit(0 if worst <= tol else 1)


if __name__ == "__main__":
    main()
