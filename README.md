# ASL

## Student : François Costa | 19-931-989


# Question 1



### Miss/hit pattern iteration 1

| Line in code | Pattern |
|---|---|
| Line 12 | MMHH |
| Line 13 | MMHH |
| Line 14 | HMHH |

### State of the cache after iteration 1

| Set | Block 0   | Block 1 |
|---|---|---|
| 1 | x(0).a x(0).b | x(3).b x(3).u(0) |
| 2 |  - | -  |
| 3 |  y(0).a y(0).b | y(3).b y(3).u(0)  |
| 4 |  x(2).u(2) x(3).a | -   |


### Miss/hit pattern iteration 2

| Line in code | Pattern |
|---|---|
| Line 12 | MMHM |
| Line 13 | MMMM |
| Line 14 | MHHM |


### State of the cache after iteration 2

| Set | Block 0   | Block 1 |
|---|---|---|
| 1 | x(0).a x(0).b | x(3).b x(3).u(0) |
| 2 |  x(2).a x(2).b | -  |
| 3 |  x(0).u(2) x(1).a | y(3).b y(3).u(0)  |
| 4 |  x(2).u(2) x(3).a | y(2).a y(2).b |

### Miss/hit pattern iteration 3

| Line in code | Pattern |
|---|---|
| Line 12 | HHHM |
| Line 13 | MHHM |
| Line 14 | HHHM |

### State of the cache after iteration 3

| Set | Block 0   | Block 1 |
|---|---|---|
| 1 | x(3).b x(3).u(0)  | x(0).a x(0).b |
| 2 |  x(2).a x(2).b | y(1).b y(1).u(0)  |
| 3 |  x(0).u(2) x(1).a | y(3).b y(3).u(0)  |
| 4 |  x(1).b x(1).u(0) | y(2).a y(2).b |


