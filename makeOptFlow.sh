# Specify the path to the optical flow utility here.
# Also check line 44 and 47 whether the arguments are in the correct order.
flowCommandLine="bash run-deepflow.sh"

if [ -z "$flowCommandLine" ]; then
  echo "Please open makeOptFlow.sh and specify the command line for computing the optical flow."
  exit 1
fi

if [ ! -f ./consistencyChecker/consistencyChecker ]; then
  if [ ! -f ./consistencyChecker/Makefile ]; then
    echo "Consistency checker makefile not found."
    exit 1
  fi
  cd consistencyChecker/
  make
  cd ..
fi

filePattern=$1
folderName=$2
startFrame=${3:-1}
stepSize=${4:-1}

if [ "$#" -le 1 ]; then
   echo "Usage: ./makeOptFlow <filePattern> <outputFolder> [<startNumber> [<stepSize>]]"
   echo -e "\tfilePattern:\tFilename pattern of the frames of the videos."
   echo -e "\toutputFolder:\tOutput folder."
   echo -e "\tstartNumber:\tThe index of the first frame. Default: 1"
   echo -e "\tstepSize:\tThe step size to create long-term flow. Default: 1"
   exit 1
fi

i=$[$startFrame]
j=$[$startFrame + $stepSize]
# j might be less than i?

mkdir -p "${folderName}"

 while true; do
  file1=$(printf "$filePattern" "$i")
  file2=$(printf "$filePattern" "$j")
  if [ -a $file2 ]; then
    if [ ! -f ${folderName}/forward_${i}_${j}.flo ]; then
      eval $flowCommandLine "$file1" "$file2" "${folderName}/forward_${i}_${j}.flo"
    fi
    if [ ! -f ${folderName}/backward_${j}_${i}.flo ]; then
      eval $flowCommandLine "$file2" "$file1" "${folderName}/backward_${j}_${i}.flo"
    fi
    # write out for both directions (backward->forward, forward->backward)
    # note j means something like "temporal prior".
    # for a backward->forward pass, this is a negative number.
    # for a forward->backward pass, this is a positive number

    # backward->forward
    # backwards_1011_1010.flo, forward_1010_1011.flo, disocclusion_1011_1010
    # prman style: backwards_1011.exr, forward_1010.exr, i.e. fore1010 and back 1011 are a pair
    # note the same pair is used for the disocclusion mask in both directions
    ./consistencyChecker/consistencyChecker "${folderName}/backward_${j}_${i}.flo" "${folderName}/forward_${i}_${j}.flo" "${folderName}/reliable_${j}_${i}.pgm"

    ./consistencyChecker/consistencyChecker "${folderName}/forward_${i}_${j}.flo" "${folderName}/backward_${j}_${i}.flo" "${folderName}/reliable_${i}_${j}.pgm"
  else
    break
  fi
  i=$[$i +1]
  j=$[$j +1]
done
