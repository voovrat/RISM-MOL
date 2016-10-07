function R = unit2unit(Tab,unit1,unit2)

unit1(unit1=='/') = '_';
unit2(unit2=='/') = '_';

V1 = eval( [ 'Tab.' unit1]);
V2 = eval( [ 'Tab.' unit2]);

R = V1/V2;