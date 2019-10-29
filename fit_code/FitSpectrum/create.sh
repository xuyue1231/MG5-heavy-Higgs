for((i=2; i<4; i++)); do
  let down=($i-1)*5+1
  let up=($i-1)*5+7
  cp fit_com_1.cc fit_com_$i.cc
  sed -i "s/alpha_1/alpha_$i/g" fit_com_$i.cc
  sed -i "s/Row>1/Row>$down/g" fit_com_$i.cc
  sed -i "s/Row<7/Row<$up/g" fit_com_$i.cc
  cp submit_1.sh submit_$i.sh
  sed -i "s/out1/out$i/g" submit_$i.sh
  cp out1.cc out$i.cc
  sed -i "s/out1/out$i/g" out$i.cc
  sed -i "s/com_1/com_$i/g" out$i.cc
done
