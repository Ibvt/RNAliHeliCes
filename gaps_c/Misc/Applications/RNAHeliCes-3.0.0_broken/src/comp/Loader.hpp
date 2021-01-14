//*****************************************************************************
#ifndef __LOADER_H
#define __LOADER_H 1

/**
@file Loader.h
**************
@brief Loader.h is the base class header for loaders.

Loader.h is the base class header for bioinformatics structures loaders.
It defines the namespace val::biocpp::loaders.

@see val::biocpp::loaders

@author Valentin GUIGNON
@version 1.0
@date 06/09/2004
*/
namespace val
{
namespace biocpp
{
/**
@namespace val::biocpp::loaders
*******************************
@brief loaders namespace contains various loaders of bioinformatics structures.

@c loaders is a namespace that contains various loaders for bioinformatics
structures such as primary structures, secondary structures. Loaders can
load a structure from various sources like files or strings.
*/
namespace loaders
{

/**
@struct SLoader   Loader.h "Loader.hpp"
***************
@brief SLoader is a template for bioinformatics structures loaders.

@c SLoader is a template for bioinformatics structures loaders. It can't be
instanciated and only describes what a loader must implement in order to
work properly. Loaders are just a kind of functors that can instanciate a
bioinformatics structure from a given data. They can be used to load structures
from files, strings, databases, caches and more.

@b Templates: \n
@param Type: \n
 the type of the bioinformatics structure to return.

@param Key: \n
 the type of the key that can identify the structure to load. For example, to
 load a structure from a file, this could be either a file type or a string
 that contains the file path and name.
*/
template <class Type, class Key>
struct SLoader
{
	virtual Type operator()(const Key &k) = 0;
	virtual ~SLoader(void) {}
}; // struct SLoader

}; // namespace loaders
}; // namespace biocpp
}; // namespace val
#endif //ifndef __LOADER_H
