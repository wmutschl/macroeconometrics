'Requires galiqje.wf1

pageselect gali_92

smpl 1955q1 1987q3

delete(noerr) gali_sys

system gali_sys
gali_sys.append YGR = C(1)*YGR(-1) + C(2)*YGR(-2) + C(3)*YGR(-3) + C(4)*YGR(-4) + C(5)*DRATE(-1) + C(6)*DRATE(-2) + C(7)*DRATE(-3) + C(8)*DRATE(-4) + C(9)*EC1(-1) + C(10)*EC1(-2) + C(11)*EC1(-3) + C(12)*EC1(-4) + C(13)*EC2(-1) + C(14)*EC2(-2) + C(15)*EC2(-3) + C(16)*EC2(-4) + C(17) - (C(5)+C(6)+C(7)+C(8))*DRATE - (C(9)+C(10)+C(11)+C(12))*EC1 - (C(13)+C(14)+C(15)+C(16))*EC2
gali_sys.append DRATE = C(18)*YGR(-1) + C(19)*YGR(-2) + C(20)*YGR(-3) + C(21)*YGR(-4) + C(22)*DRATE(-1) + C(23)*DRATE(-2) + C(24)*DRATE(-3) + C(25)*DRATE(-4) + C(26)*EC1(-1) + C(27)*EC1(-2) + C(28)*EC1(-3) + C(29)*EC1(-4) + C(30)*EC2(-1) + C(31)*EC2(-2) + C(32)*EC2(-3) + C(33)*EC2(-4) + C(34) + C(69)*YGR + C(70)*(EC1-EC2)
gali_sys.append EC1 = C(35)*YGR(-1) + C(36)*YGR(-2) + C(37)*YGR(-3) + C(38)*YGR(-4) + C(39)*DRATE(-1) + C(40)*DRATE(-2) + C(41)*DRATE(-3) + C(42)*DRATE(-4) + C(43)*EC1(-1) + C(44)*EC1(-2) + C(45)*EC1(-3) + C(46)*EC1(-4) + C(47)*EC2(-1) + C(48)*EC2(-2) + C(49)*EC2(-3) + C(50)*EC2(-4) + C(51) + C(71)*YGR + C(72)*DRATE + C(73)*EC2
gali_sys.append EC2 = C(52)*YGR(-1) + C(53)*YGR(-2) + C(54)*YGR(-3) + C(55)*YGR(-4) + C(56)*DRATE(-1) + C(57)*DRATE(-2) + C(58)*DRATE(-3) + C(59)*DRATE(-4) + C(60)*EC1(-1) + C(61)*EC1(-2) + C(62)*EC1(-3) + C(63)*EC1(-4) + C(64)*EC2(-1) + C(65)*EC2(-2) + C(66)*EC2(-3) + C(67)*EC2(-4) + C(68) + C(74)*YGR -((C(5)+C(6)+C(7)+C(8))/(C(13)+C(14)+C(15)+C(16)))*DRATE -((C(9)+C(10)+C(11)+C(12))/(C(13)+C(14)+C(15)+C(16)))*EC1
gali_sys.ls
gali_sys.results


