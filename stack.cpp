#include <iostream>
#include "stack.h"

using namespace std;

Stack::Stack()
{
    up = NULL;
}

void Stack::push(std::string nameParam)
{
    Cell *newCell = new Cell;

    newCell->after = up;
    newCell->name = nameParam;

    up = newCell;
}

std::string Stack::peek()
{
    if(!isEmpty())
    {
        return up->name;
    }

    return 0;
}

std::string Stack::pop()
{
    if(!isEmpty())
    {
        std::string name;
        Cell *cell = NULL;

        name = up->name;

        cell = up->after;

        delete up;

        up = cell;

        return name;
    }

    return 0;
}

bool Stack::isEmpty() const
{
    if(up == NULL)
    {
        return true;
    }
    return false;
}

void Stack::clear()
{
    while(!isEmpty())
    {
        pop();
    }
}

Stack::~Stack()
{
    clear();
}
