# Structuring a python project

In this tutorial you will learn how to structure a python project.

We'll create a small python package step by step.
Our package will convert mass to weight and vice-versa.

!!! note
    `nano` will start to show its limits as a file editor. We recommend the excellent `atom` for this lesson.
    Additionally it has built-in git support!

!!! warning
    Don't forget to commit your changes on a regular basis

## Creating your project

We will call you project `newton`

```bash
mkdir newton
cd newton
```

## The structure of a project

Inside your newly created `newton` directory, our finished python project will look like this:

```
.
..
newton/
    __init__.py
    __main__.py
    app.py
    converter.py
    test/
        test_converter.py
setup.py
LICENSE
README.md
```

Don't worry, we'll explain what each file does, and create them along the way

## Modules

Our `newton` package will be both a library (you will be able to do `import newton`) and a command-line application (you will be able to type `newton` in your terminal)

For that we need a few things. Let us create a directory called `newton` within our already existing newtonon directory

### \_\_init\_\_.py

```python
from newton import *
```

### \_\_main\_\_.py

```python
from newton.app import main
main()
```

### app.py

This is the core of our project. it should contain a `main` function

```python
def main():
    print('Hello! My name is newton!')
```

and now in the terminal:

```bash
python -m newton
# Hello! My name is newton!
```

## Command-line arguments

Python has a built-in library called `argparse` that we can use to add arguments to our package.

```python
import argparse


def main():
    parser = argparse.ArgumentParser(
        prog='newton',
        description='mass to weight converter'
    )
    parser.add_argument(
        '--mass',
        '-m',
        help='mass in kg',
        required=True
    )
    args = parser.parse_args()

    print('mass is %s kg' % args.mass)
```

we now have a help message!

```bash
python -m newton -h
```

and if we execute

```bash
python -m newton --mass 61
```

## Modular design

Every new functionality of our application should be its own function, or set of functions in a module as it grows

Let us write our mass to weight converter in a new function

```python
def mass_to_weight(mass):
    g = 9.81
    weight = mass * g
    return weight
```

and at the end of our main:

```python
print('mass is %s Kg' % args.mass)
weight = mass_to_weight(args.mass)
print('weight is %s N' % weight)
```

and if we execute

```bash
python -m newton --mass 61
```

```python-traceback
mass is 61 Kg
Traceback (most recent call last):
  File "/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/runpy.py", line 162, in _run_module_as_main
    "__main__", fname, loader, pkg_name)
  File "/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/runpy.py", line 72, in _run_code
    exec code in run_globals
  File "/Users/hadrien/Documents/workspace/newton/newton/__main__.py", line 2, in <module>
    main()
  File "newton/app.py", line 24, in main
    weight = mass_to_weight(args.mass)
  File "newton/app.py", line 6, in mass_to_weight
    weight = mass * g
TypeError: can't multiply sequence by non-int of type 'float'
```

it fails!
Indeed, the default type for argparse arguments in strings.
We need to explicitely tell argparse that we are expecting a number.

```python
parser.add_argument(
    '--mass',
    '-m',
    help='mass in Kg',
    type=float,
    required=True
)
```

now we execute again

```bash
python -m newton --mass 61
# mass is 61.0 Kg
# weight is 598.41 N
python -m newton --mass 10
# mass is 10.0 Kg
# weight is 98.10000000000001 N
```

Computers are not very good a float operations

let us modify our function to round results at two decimals

### Do not clutter `app.py`

As your application grows, more and more things will end up in your app.py.
It is good practise to move the core functionalities of your application in separate files.
It will help future you with testing, debugging and reading your code!

We create a new file, `converter.py` and move the mass_to_weight function to it.

Now, for the `mass_to_weight` function to be available in your app, you need to import it

```python
from newton import converter
```

and the line calculating the weight becomes

```python
weight = converter.mass_to_weight(args.mass)
```

## Writing tests

Now, while your function is relatively simple, we'd like to be sure it performs as expected.
Enter unit tests.

Tests are a way of making sure that any change we may make to the code doesn't break the intended functionality

### nosetests

nose is a unit test frawork for python

install it from your terminal

```bash
pip install nose
```

then we create a `test` directory in our `newton` module directory

```bash
mkdir -p newton/test
```

In that test directory, we create a file named `test_converter.py`, containing

```python
from newton.converter import mass_to_weight


def test_mass_to_weight():
    weight = mass_to_weight(10)
    assert weight == 98.1
```

then in our terminal

```bash
nosetests
```

## Defensive programming

Right now, your function seem to work as intended.
Thanks to argparse, if we try to break our program by giving a string instead of a float:

```bash
python -m newton --mass twenty
# usage: newton [-h] --mass MASS
# newton: error: argument --mass/-m: invalid float value: 'twenty'
```

But what happens if we enter a negative value?

```bash
python -m newton --mass -12
# mass is -12.0 Kg
# weight is -117.72 N
```

That doesn't seem right.
We should restrict our function to only run with positive mass and display an error when the mass is negative

```python
def mass_to_weight(mass):
    g = 9.81
    try:
        assert mass >= 0
        weight = mass * g
    except AssertionError as e:
        raise e
    else:
        weight = mass * g
        return round(weight, 2)
```

If we run our app again with a negative mass, we now have an AssertionError, accompanied by a traceback.
While effective, this is not very user friendly

```python
import sys


def mass_to_weight(mass):
    g = 9.81
    try:
        assert mass >= 0
        weight = mass * g
    except AssertionError as e:
        print('mass cannot be negative!')
        sys.exit(1)
    else:
        weight = mass * g
        return round(weight, 2)
```

The code snippet above on the other hand, will print a message explaining the error, and then exit with the code 1.

!!! note
     This is a convention for command line programs to exit with 0 when they succeed, and 1 or another error code when they fail.

If we run our test again, it should still pass.
Let us write another test that makes sure we handle negative values correctly

In test_converter.py

```python
from nose.tools import raises

# ...

@raises(SystemExit)
def test_negative_mass():
    weight = mass_to_weight(-10)
```

```bash
nosetests
```

should show that the two tests have passed.

!!! question
    What happens when the mass is zero?
    Write a test to make sure our code handles that case properly

## Adding functionalities

Now our python program does only one thing, even it it does it well.
We'll add the following functionalities to our program.

* Earth is not the only planet that matters! It should be possible to change the g parameter (but still have 9.81 as a default value)

* It should be possible to convert weight to mass

## The README

## A word on Licenses
