# Time-series benchmark data

These fixed CSV files are real-world univariate series from R's base
`datasets` package, converted to CSV by the Rdatasets archive. They are kept
small enough to live in the repository so benchmark runs require no network
access.

The files are pinned to Rdatasets commit
`1dcc2bf5f955cc1224a3e1307256e1fe86b68dae`. To replace a file reproducibly,
download:

```text
https://raw.githubusercontent.com/vincentarelbundock/Rdatasets/1dcc2bf5f955cc1224a3e1307256e1fe86b68dae/csv/datasets/<name>.csv
```

| File | Values | Period | Benchmark role | SHA-256 |
| --- | ---: | ---: | --- | --- |
| `AirPassengers.csv` | 144 | 12 | Integrated seasonal airline demand | `bdb98adbd418a6de6842a742e0602f363c3b841c26677e7e61a1e055e9509bd8` |
| `Nile.csv` | 100 | 1 | Short AR and exact-likelihood fit | `d0452bea38c61e796a4eeb950bf91d20fb5c7f13d5822eadf5990fe54f9c8d07` |
| `nottem.csv` | 240 | 12 | Stationary seasonal temperature | `634e2e374c5d59c1332e6e6de207d691cb4ebdcd007f7f0e2993ac2a19c90473` |
| `sunspots.csv` | 2,820 | 12 | Long ARMA fit | `8db04db1d4406f8a005c00fc3a4647de1845fd3e7e7bbde4ad176a19037c9e60` |
| `UKgas.csv` | 108 | 4 | Integrated quarterly seasonality | `69ec013fbeeaf956a42603aaf8f68d943a3aab5f90b6ed89c5de2ed02eba59ff` |
| `WWWusage.csv` | 100 | 1 | Short integrated series and order search | `245774828aa5bc64292f29db15434004b7472022b73b81a3d5d1c02f55c42e08` |

The first column is a one-based row number, the second is R's numeric time,
and the third is the observed value. All files have the header
`rownames,time,value`.

## Provenance and licensing

Rdatasets distributes its repository code and copied documentation under
GPL-3.0. Its maintainer notes that the licensing of some underlying numeric
datasets is not definitive, while documenting a good-faith understanding that
they are redistributable. JLinAlg therefore preserves the original names and
provenance rather than asserting a new license over the observations. See the
[Rdatasets license](https://github.com/vincentarelbundock/Rdatasets/blob/1dcc2bf5f955cc1224a3e1307256e1fe86b68dae/LICENSE.md)
and [licensing note](https://github.com/vincentarelbundock/Rdatasets/tree/1dcc2bf5f955cc1224a3e1307256e1fe86b68dae#license).

Original R documentation and source citations:

- [AirPassengers](https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/AirPassengers.html)
- [Nile](https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/Nile.html)
- [nottem](https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/nottem.html)
- [sunspots](https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/sunspots.html)
- [UKgas](https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/UKgas.html)
- [WWWusage](https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/WWWusage.html)
