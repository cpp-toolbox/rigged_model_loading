#include "rigged_model_loading.hpp"
#include <assimp/postprocess.h>
#include <glm/gtx/matrix_decompose.hpp>
#include <glm/gtx/quaternion.hpp> // For quaternion operations like eulerAngles
#include <regex>
#include <filesystem>
/**
 * @brief Extracts the scale vector from a 4x4 transformation matrix.
 *
 * @param transformation_matrix The 4x4 transformation matrix.
 * @return glm::vec3 The scale factors along the X, Y, and Z axes.
 */
glm::vec3 parse_scale_matrix(const glm::mat4 &transformation_matrix) {
    glm::vec3 scale;
    scale.x = glm::length(glm::vec3(transformation_matrix[0])); // Length of the X-axis column vector
    scale.y = glm::length(glm::vec3(transformation_matrix[1])); // Length of the Y-axis column vector
    scale.z = glm::length(glm::vec3(transformation_matrix[2])); // Length of the Z-axis column vector
    return scale;
}

/**
 * @brief Prints the details of an aiAnimation instance.
 *
 * @param animation The aiAnimation instance to print.
 */
void print_ai_animation_info(const aiAnimation *animation) {
    if (!animation) {
        std::cout << "Invalid animation pointer.\n";
        return;
    }

    // Print basic information
    std::cout << "Animation Name: " << animation->mName.C_Str() << "\n";
    std::cout << "Duration: " << animation->mDuration << " ticks\n";
    std::cout << "Ticks per second: "
              << (animation->mTicksPerSecond == 0 ? "Not specified" : std::to_string(animation->mTicksPerSecond))
              << "\n";

    // Print node animation channels
    std::cout << "Number of Node Animation Channels: " << animation->mNumChannels << "\n";
    for (unsigned int i = 0; i < animation->mNumChannels; ++i) {
        const aiNodeAnim *channel = animation->mChannels[i];
        if (channel) {
            std::cout << "  Channel " << i + 1 << " Node Name: " << channel->mNodeName.C_Str() << "\n";
        }
    }

    // Print mesh animation channels
    std::cout << "Number of Mesh Animation Channels: " << animation->mNumMeshChannels << "\n";
    for (unsigned int i = 0; i < animation->mNumMeshChannels; ++i) {
        const aiMeshAnim *meshChannel = animation->mMeshChannels[i];
        if (meshChannel) {
            std::cout << "  Mesh Channel " << i + 1 << " Mesh Name: " << meshChannel->mName.C_Str() << "\n";
        }
    }

    // Print morph mesh animation channels
    std::cout << "Number of Morph Mesh Animation Channels: " << animation->mNumMorphMeshChannels << "\n";
    for (unsigned int i = 0; i < animation->mNumMorphMeshChannels; ++i) {
        const aiMeshMorphAnim *morphChannel = animation->mMorphMeshChannels[i];
        if (morphChannel) {
            std::cout << "  Morph Channel " << i + 1 << " Mesh Name: " << morphChannel->mName.C_Str() << "\n";
        }
    }
}

/**
 * @brief Extracts the translation vector from a 4x4 transformation matrix.
 *
 * @param transformation_matrix The 4x4 transformation matrix.
 * @return glm::vec3 The translation vector (X, Y, Z).
 */
glm::vec3 parse_translation_matrix(const glm::mat4 &transformation_matrix) {
    return glm::vec3(transformation_matrix[3][0], transformation_matrix[3][1], transformation_matrix[3][2]);
}

/**
 * @brief Converts a rotation matrix to Euler angles (in degrees).
 *
 * @param rotation_matrix The 4x4 rotation matrix.
 * @return glm::vec3 A vector containing Euler angles (pitch, yaw, roll) in degrees.
 */
glm::vec3 matrix_to_euler_angles(const glm::mat4 &rotation_matrix) {
    // Extract the rotation part (upper-left 3x3 matrix)
    glm::mat3 rot(rotation_matrix);

    // Calculate Euler angles
    float sy = std::sqrt(rot[0][0] * rot[0][0] + rot[1][0] * rot[1][0]);

    bool singular = sy < 1e-6; // Check for gimbal lock (singular case)

    float pitch, yaw, roll;
    if (!singular) {
        pitch = std::atan2(rot[2][1], rot[2][2]);
        yaw = std::atan2(-rot[2][0], sy);
        roll = std::atan2(rot[1][0], rot[0][0]);
    } else {
        pitch = std::atan2(-rot[1][2], rot[1][1]);
        yaw = std::atan2(-rot[2][0], sy);
        roll = 0; // Roll is undefined in singularity; we set it to 0
    }

    // Convert from radians to degrees
    return glm::degrees(glm::vec3(pitch, yaw, roll));
}

/**
 * @brief Print a vector (scale, translation, or Euler angles) in a readable format.
 *
 * @param vec The glm::vec3 containing the values to print.
 * @param description A description of what the vector represents.
 */
void print_vector(const glm::vec3 &vec, const std::string &description, int indentation_level = 0) {
    std::string indentation(indentation_level * 2, '   '); // 2 spaces per level
    std::cout << indentation << description << ":" << std::endl;
    std::cout << indentation << "X: " << vec.x << std::endl;
    std::cout << indentation << "Y: " << vec.y << std::endl;
    std::cout << indentation << "Z: " << vec.z << std::endl;
}

void print_matrix(const glm::mat4 &matrix, const std::string &description, int indentation_level = 0) {
    // Generate indentation based on the provided level
    std::string indentation(indentation_level * 2, '   '); // 2 spaces per level

    // Print description with indentation
    std::cout << indentation << description << ":\n";

    // Print top border of the box
    std::cout << indentation << "+------------------------+" << std::endl;

    // Print matrix rows with box around them
    for (int row = 0; row < 4; ++row) {
        std::cout << indentation << "| "; // Left border
        for (int col = 0; col < 4; ++col) {
            std::cout << matrix[row][col] << " "; // Matrix values
        }
        std::cout << std::endl; // Right border
    }

    // Print bottom border of the box
    std::cout << indentation << "+------------------------+" << std::endl;
}

void print_euler_angles(const glm::quat &quat, const std::string &description, int indentation_level = 0) {
    std::string indentation(indentation_level * 2, ' ');           // 2 spaces per level
    glm::vec3 euler_angles = glm::degrees(glm::eulerAngles(quat)); // Convert to degrees
    std::cout << indentation << description << " (Euler Angles):" << std::endl;
    std::cout << indentation << "  Pitch (X): " << euler_angles.x << "°" << std::endl;
    std::cout << indentation << "  Yaw (Y): " << euler_angles.y << "°" << std::endl;
    std::cout << indentation << "  Roll (Z): " << euler_angles.z << "°" << std::endl;
}

/**
 * @brief Decomposes and prints a 4x4 transform matrix with a description.
 *
 * @param matrix The 4x4 transformation matrix (glm::mat4).
 * @param description A description of the matrix (e.g., "Transform Matrix").
 * @param indentation_level Indentation level for nested structures.
 */
void print_transform(const glm::mat4 &matrix, const std::string &description, int indentation_level = 0) {
    std::string indentation(indentation_level * 2, ' '); // 2 spaces per level

    std::cout << indentation << "+------START-------+" << std::endl;
    std::cout << indentation << description << ":" << std::endl;

    glm::vec3 scale, translation, skew;
    glm::quat rotation;
    glm::vec4 perspective;

    // Use glm::decompose to extract transform components
    glm::decompose(matrix, scale, rotation, translation, skew, perspective);

    // Print components
    print_vector(translation, "Translation", indentation_level + 1);
    print_vector(scale, "Scale", indentation_level + 1);
    print_euler_angles(rotation, "Rotation", indentation_level + 1);

    std::cout << indentation << "+-------END--------+" << std::endl;
}

glm::mat4 ai_matrix4x4_to_glm_mat4(const aiMatrix4x4 &ai_mat) {
    glm::mat4 glm_mat;

    // Transpose the row-major aiMatrix4x4 to column-major glm::mat4
    glm_mat[0][0] = ai_mat.a1;
    glm_mat[1][0] = ai_mat.a2;
    glm_mat[2][0] = ai_mat.a3;
    glm_mat[3][0] = ai_mat.a4;
    glm_mat[0][1] = ai_mat.b1;
    glm_mat[1][1] = ai_mat.b2;
    glm_mat[2][1] = ai_mat.b3;
    glm_mat[3][1] = ai_mat.b4;
    glm_mat[0][2] = ai_mat.c1;
    glm_mat[1][2] = ai_mat.c2;
    glm_mat[2][2] = ai_mat.c3;
    glm_mat[3][2] = ai_mat.c4;
    glm_mat[0][3] = ai_mat.d1;
    glm_mat[1][3] = ai_mat.d2;
    glm_mat[2][3] = ai_mat.d3;
    glm_mat[3][3] = ai_mat.d4;

    return glm_mat;
}

glm::mat4 ai_matrix3x3_to_glm_mat4(const aiMatrix3x3 &ai_mat) {
    // Create a glm::mat4 matrix initialized to the identity matrix
    glm::mat4 glmMat(1.0f);

    // Set the upper-left 3x3 part of the glm::mat4 matrix
    glmMat[0][0] = ai_mat.a1;
    glmMat[1][0] = ai_mat.a2;
    glmMat[2][0] = ai_mat.a3;
    glmMat[0][1] = ai_mat.b1;
    glmMat[1][1] = ai_mat.b2;
    glmMat[2][1] = ai_mat.b3;
    glmMat[0][2] = ai_mat.c1;
    glmMat[1][2] = ai_mat.c2;
    glmMat[2][2] = ai_mat.c3;

    return glmMat;
}

void print_ai_animation(const aiAnimation *anim) {
    // Draw top border
    std::cout << "+-----------------------------------------------+" << std::endl;

    std::cout << "| Animation Name: " << anim->mName.data << std::endl;
    std::cout << "| Duration: " << anim->mDuration << " ticks" << std::endl;
    std::cout << "| Ticks per second: " << anim->mTicksPerSecond << std::endl;
    std::cout << "| Number of Bone Animation Channels: " << anim->mNumChannels << std::endl;
    std::cout << "| Number of Mesh Animation Channels: " << anim->mNumMeshChannels << std::endl;
    std::cout << "| Number of Morph Mesh Animation Channels: " << anim->mNumMorphMeshChannels << std::endl;

    // Printing Bone Animation Channels
    if (anim->mChannels) {
        std::cout << "| Bone Animation Channels: " << std::endl;
        for (unsigned int i = 0; i < anim->mNumChannels; ++i) {
            std::cout << "|   Channel " << i + 1 << ": " << &(anim->mChannels[i]) << std::endl;
        }
    } else {
        std::cout << "| No Bone Animation Channels" << std::endl;
    }

    // Printing Mesh Animation Channels
    if (anim->mMeshChannels) {
        std::cout << "| Mesh Animation Channels: " << std::endl;
        for (unsigned int i = 0; i < anim->mNumMeshChannels; ++i) {
            std::cout << "|   Channel " << i + 1 << ": " << &(anim->mMeshChannels[i]) << std::endl;
        }
    } else {
        std::cout << "| No Mesh Animation Channels" << std::endl;
    }

    // Printing Morph Mesh Animation Channels
    if (anim->mMorphMeshChannels) {
        std::cout << "| Morph Mesh Animation Channels: " << std::endl;
        for (unsigned int i = 0; i < anim->mNumMorphMeshChannels; ++i) {
            std::cout << "|   Channel " << i + 1 << ": " << &(anim->mMorphMeshChannels[i]) << std::endl;
        }
    } else {
        std::cout << "| No Morph Mesh Animation Channels" << std::endl;
    }

    // Draw bottom border
    std::cout << "+-----------------------------------------------+" << std::endl;
}

/**
 * @brief Prints the mapping of armature names to animation indices.
 *
 * @param map A map of armature names to their corresponding animation indices.
 */
void print_armature_to_animation_map(const std::unordered_map<std::string, unsigned int> &map) {
    std::cout << "Armature to Animation Mapping:" << std::endl;

    for (const auto &pair : map) {
        const std::string &armature_name = pair.first;
        unsigned int animation_index = pair.second;

        std::cout << "Armature: " << armature_name << std::endl;
        std::cout << "  Animation Index: " << animation_index << std::endl;
    }
}

// should only ever get called once
void RecIvpntRiggedCollector::rec_process_nodes(aiNode *node, const aiScene *scene) {
    // Helper to generate indentation based on the recursion level
    auto get_indentation = [this]() {
        return std::string(recursion_level_counter * 2, ' '); // 2 spaces per level
    };

    bool logging = false;

    if (logging) {
        std::cout << get_indentation() << "Recursion Level: " << recursion_level_counter << std::endl;

        // Print the name of the current node
        if (node->mName.length > 0) {
            std::cout << get_indentation() << "Processing Node: " << node->mName.C_Str() << std::endl;
        } else {
            std::cout << get_indentation() << "Processing Node: (Unnamed)" << std::endl;
        }

        // Print the number of meshes in the current node
        std::cout << get_indentation() << "Number of Meshes in Node: " << node->mNumMeshes << std::endl;
    }

    // Process each mesh in the current node
    for (unsigned int i = 0; i < node->mNumMeshes; i++) {
        unsigned int mesh_index = node->mMeshes[i];
        aiMesh *mesh = scene->mMeshes[mesh_index];

        if (logging) {
            // Print details about the mesh
            std::cout << get_indentation() << "  Mesh Index: " << mesh_index << std::endl;
            std::cout << get_indentation()
                      << "  Mesh Name: " << (mesh->mName.length > 0 ? mesh->mName.C_Str() : "(Unnamed)") << std::endl;
            std::cout << get_indentation() << "  Number of Vertices: " << mesh->mNumVertices << std::endl;
            std::cout << get_indentation() << "  Number of Faces: " << mesh->mNumFaces << std::endl;
        }

        // Store processed mesh data
        this->ivpntrs.push_back(process_mesh_ivpntrs(mesh, scene));
    }

    if (logging) {
        // Print the number of children for the current node
        std::cout << get_indentation() << "Number of Children: " << node->mNumChildren << std::endl;
    }

    // Recurse into child nodes
    for (unsigned int i = 0; i < node->mNumChildren; i++) {
        if (logging) {
            std::cout << get_indentation() << "Recursing into Child " << i + 1 << "/" << node->mNumChildren
                      << std::endl;
        }
        recursion_level_counter++;
        this->rec_process_nodes(node->mChildren[i], scene);
        recursion_level_counter--; // Decrement after returning from recursion
    }

    if (logging) {
        // Indicate end of node processing
        std::cout << get_indentation()
                  << "Finished Processing Node: " << (node->mName.length > 0 ? node->mName.C_Str() : "(Unnamed)")
                  << std::endl;
    }
}

void RecIvpntRiggedCollector::rec_update_animation_matrices(float animation_time_ticks,
                                                            glm::mat4 parent_animation_transform_in_local_space,
                                                            aiNode *node, const aiScene *scene, int rec_depth) {
    // Helper to generate indentation based on the recursion level
    auto get_indentation = [this, &rec_depth]() {
        return std::string(rec_depth * 2, ' '); // 2 spaces per level
    };

    bool logging = false;

    if (logging) {
        std::cout << "rec_update_animation_matrices just called with animation time: " << animation_time_ticks
                  << " rec_depth: " << rec_depth << std::endl;

        print_matrix(parent_animation_transform_in_local_space, "parent_transform", rec_depth);
    }

    glm::mat4 animation_transform_for_current_time_in_bone_space = ai_matrix4x4_to_glm_mat4(node->mTransformation);

    std::string node_name(node->mName.data);

    if (logging) {
        // Indented print output
        std::cout << get_indentation() << "on node: " << node_name << std::endl;
    }

    // the following works because an armature node is the root of all the bones.
    bool node_is_armature = is_armature_node(node);
    if (node_is_armature) {
        /*curr_armature_name_rec = node_name;*/
        curr_animation_index_rec = armature_node_name_to_animation_index.at(node_name);
        if (logging) {
            // Indented print output
            std::cout << "this node is an armature using the following animation index:" << curr_animation_index_rec
                      << std::endl;
        }
    }

    // TODO in the future load different animations
    const aiAnimation *animation = scene->mAnimations[curr_animation_index_rec];
    const aiNodeAnim *node_anim = find_node_anim(animation, node_name);

    bool node_is_animated = node_anim != NULL;
    bool user_requested_no_anim = animation_time_ticks == -1;

    if (not user_requested_no_anim and node_is_animated) {

        if (logging) {
            std::cout << get_indentation() << "current node has an animation" << std::endl;
        }

        aiVector3D scaling;
        calc_interpolated_scaling(scaling, animation_time_ticks, node_anim);
        glm::vec3 glm_scaling(scaling.x, scaling.y, scaling.z);
        glm::mat4 scale_transform = glm::scale(glm::mat4(1.0f), glm_scaling);
        /*scale_transform = glm::mat4(1.0);*/

        // Indented print matrix output
        /*print_matrix(scale_transform, "scale matrix", rec_depth);*/
        /*print_vector(parse_scale_matrix(scale_transform), "scale_matrix", rec_depth);*/

        aiQuaternion rotation;
        calc_interpolated_rotation(rotation, animation_time_ticks, node_anim);
        glm::mat4 rotation_transform = ai_matrix3x3_to_glm_mat4(rotation.GetMatrix());
        // Indented print matrix output
        /*print_vector(matrix_to_euler_angles(rotation_transform), "rotation in euler angles", rec_depth);*/

        aiVector3D translation;
        calc_interpolated_translation(translation, animation_time_ticks, node_anim);
        glm::vec3 glm_translation(translation.x, translation.y, translation.z);
        glm::mat4 translation_transform = glm::translate(glm::mat4(1.0f), glm_translation);
        // Indented print matrix output
        /*print_vector(parse_translation_matrix(translation_transform), "translation", rec_depth);*/

        // we overwrite here based on assimp's documentation, when there is animation we don't use
        // mTransformation
        animation_transform_for_current_time_in_bone_space =
            translation_transform * rotation_transform * scale_transform;
        if (logging) {
            print_transform(animation_transform_for_current_time_in_bone_space, "animation transform (trs)", rec_depth);
        }
    }

    if (logging) {
        std::cout << get_indentation() << "computed matrices" << std::endl;
    }

    // note that the recursion goes outward towards leaves, but we think the other way, associativity of
    // matrix multiplication reconciles this.
    glm::mat4 bone_to_local_animation_transform_up_to_this_node =
        parent_animation_transform_in_local_space * animation_transform_for_current_time_in_bone_space;

    if (logging) {
        print_matrix(parent_animation_transform_in_local_space, "parent_transform", rec_depth);
        print_matrix(animation_transform_for_current_time_in_bone_space, " animation_transform_for_current_time",
                     rec_depth);
        print_transform(animation_transform_for_current_time_in_bone_space, " animation_transform_for_current_time");
        print_matrix(bone_to_local_animation_transform_up_to_this_node, "bone_to_local_transform_up_to_this_node ",
                     rec_depth);
    }

    // this variable describes the requirement for membership into the bone_name_to_unique_index
    // the vector is constructed when parse_model_into_ivptrs is called
    /*bool this_node_is_a_bone_and_weve_processed_its_vertex_data =*/
    bool node_is_bone = bone_name_to_unique_index.find(node_name) != bone_name_to_unique_index.end();

    if (node_is_bone) {

        if (logging) {
            std::cout << get_indentation() << "this node was a bone" << std::endl;
        }
        int bone_idx = bone_name_to_unique_index[node_name];
        auto &bi =
            bone_unique_idx_to_info[bone_idx]; // this is guaranteed safe cause already exists in there for some reason

        // the idea here is to first go to the bone coordinate system, then apply all the transformations that should
        // be applied due to the recursion up to this bone, and also make sure that if the armature is displaced, it
        // also moves via the inverse root_note_transform
        bi.local_space_animated_transform_upto_this_bone = inverse_root_node_transform *
                                                           bone_to_local_animation_transform_up_to_this_node *
                                                           bi.local_space_to_bone_space_in_bind_pose_transformation;

        if (logging) {
            print_transform(inverse_root_node_transform, "inverse_root_node_transform", rec_depth);
            print_transform(bone_to_local_animation_transform_up_to_this_node,
                            "bone_to_local_transform_up_to_this_node", rec_depth);
            print_transform(bi.local_space_to_bone_space_in_bind_pose_transformation,
                            "bi.local_space_to_bone_space_in_bind_pose_transformation", rec_depth);

            print_transform(bi.local_space_animated_transform_upto_this_bone,
                            "full_bone_space_to_local_space_transformation", rec_depth);
        }
    } else {
        if (logging) {
            std::cout << get_indentation() << "this node was NOT a bone" << std::endl;
        }
    }

    glm::mat4 curr_mat;
    if (node_is_bone) {
        curr_mat = bone_to_local_animation_transform_up_to_this_node;
    } else {
        // pass it through if not a bone could be bad
        /*curr_mat = parent_transform;*/
        curr_mat = bone_to_local_animation_transform_up_to_this_node;
    }

    /*spdlog::get(Systems::asset_loading)->info("finished processing meshes");*/
    for (unsigned int i = 0; i < node->mNumChildren; i++) {
        rec_update_animation_matrices(animation_time_ticks, curr_mat, node->mChildren[i], scene, rec_depth + 1);
    }
}

/*
 * @pre the asset of interest has been loaded already.
 */
void RecIvpntRiggedCollector::set_bone_transforms(float time_in_seconds, std::vector<glm::mat4> &transforms_to_be_set) {
    transforms_to_be_set.resize(bone_unique_idx_to_info.size());

    /*print_ai_animation(scene->mAnimations[0]);*/
    // uses 25 fps if ticks per second was not specified
    float ticks_per_second =
        (float)(scene->mAnimations[0]->mTicksPerSecond != 0 ? scene->mAnimations[0]->mTicksPerSecond : 25.0f);
    float time_in_ticks = ticks_per_second * time_in_seconds;
    float animation_time_ticks = fmod(time_in_ticks, (float)scene->mAnimations[0]->mDuration);

    bool logging = false;

    if (logging) {
        std::cout << "=== STARTING UPDATE ANIMATION MATRICES ===" << std::endl;
    }
    update_animation_matrices(animation_time_ticks);
    if (logging) {
        std::cout << "=== ENDING UPDATE ANIMATION MATRICES ===" << std::endl;
    }

    /*spdlog::info("bone info size", bone_info.size());*/
    for (unsigned int i = 0; i < bone_unique_idx_to_info.size(); i++) {
        /*spdlog::info("setting transform {}", bone_info[i].full_bone_space_to_local_space_transformation[0][0]);*/
        transforms_to_be_set[i] = bone_unique_idx_to_info[i].local_space_animated_transform_upto_this_bone;
    }
}

void calc_interpolated_scaling(aiVector3D &out, float animation_time_ticks, const aiNodeAnim *node_anim) {
    // we need at least two values to interpolate...
    if (node_anim->mNumScalingKeys == 1) {
        std::cout << "there is only one scaling key, scaling animation will not be applied" << std::endl;
        out = node_anim->mScalingKeys[0].mValue;
        return;
    }

    unsigned int scaling_idx = find_idx_of_scaling_key_for_given_time(animation_time_ticks, node_anim);
    unsigned int next_scaling_idx = scaling_idx + 1;

    assert(next_scaling_idx < node_anim->mNumScalingKeys);

    float t1 = (float)node_anim->mScalingKeys[scaling_idx].mTime;
    float t2 = (float)node_anim->mScalingKeys[next_scaling_idx].mTime;

    float delta_time = t2 - t1;

    // t1 < a_t_t < t2, so this is non-negative and works correclty
    float factor = (animation_time_ticks - (float)t1) / delta_time;
    assert(factor >= 0.0f && factor <= 1.0f);

    const aiVector3D &start_scale = node_anim->mScalingKeys[scaling_idx].mValue;
    const aiVector3D &end_scale = node_anim->mScalingKeys[next_scaling_idx].mValue;
    aiVector3D scale_delta = end_scale - start_scale;
    out = start_scale + factor * scale_delta;
}

void calc_interpolated_rotation(aiQuaternion &out, float animation_time_ticks, const aiNodeAnim *node_anim) {
    // we need at least two values to interpolate...
    if (node_anim->mNumRotationKeys == 1) {
        std::cout << "there is only one rotation key, scaling animation will not be applied" << std::endl;
        out = node_anim->mRotationKeys[0].mValue;
        return;
    }

    unsigned int rotation_idx = find_idx_of_rotation_key_for_given_time(animation_time_ticks, node_anim);
    unsigned int next_rotation_idx = rotation_idx + 1;
    assert(next_rotation_idx < node_anim->mNumRotationKeys);

    float t1 = (float)node_anim->mRotationKeys[rotation_idx].mTime;
    float t2 = (float)node_anim->mRotationKeys[next_rotation_idx].mTime;

    float delta_time = t2 - t1;
    float factor = (animation_time_ticks - t1) / delta_time;
    assert(factor >= 0.0f && factor <= 1.0f);

    const aiQuaternion &start_rotation = node_anim->mRotationKeys[rotation_idx].mValue;
    const aiQuaternion &end_rotation = node_anim->mRotationKeys[next_rotation_idx].mValue;

    aiQuaternion::Interpolate(out, start_rotation, end_rotation, factor);
    out.Normalize();
}

void calc_interpolated_translation(aiVector3D &out, float animation_time_ticks, const aiNodeAnim *node_anim) {
    // we need at least two values to interpolate...
    if (node_anim->mNumPositionKeys == 1) {
        std::cout << "there is only one position key, scaling animation will not be applied" << std::endl;
        out = node_anim->mPositionKeys[0].mValue;
        return;
    }

    unsigned int translation_idx = find_idx_of_translation_key_for_given_time(animation_time_ticks, node_anim);
    unsigned int next_translation_idx = translation_idx + 1;
    assert(next_translation_idx < node_anim->mNumPositionKeys);

    float t1 = (float)node_anim->mPositionKeys[translation_idx].mTime;
    float t2 = (float)node_anim->mPositionKeys[next_translation_idx].mTime;

    float delta_time = t2 - t1;
    float factor = (animation_time_ticks - t1) / delta_time;
    assert(factor >= 0.0f && factor <= 1.0f);

    const aiVector3D &start_translation = node_anim->mPositionKeys[translation_idx].mValue;
    const aiVector3D &end_translation = node_anim->mPositionKeys[next_translation_idx].mValue;

    aiVector3D translation_delta = end_translation - start_translation;
    out = start_translation + factor * translation_delta;
}

// Helper function to build the full path of a node
std::string get_full_node_path(const aiNode *node) {
    if (!node->mParent) {
        return node->mName.C_Str(); // Root node has no parent
    }
    return get_full_node_path(node->mParent) + "/" + node->mName.C_Str();
}

/**
 * @brief Builds a mapping of armature names to animation indices based on animation names in the given scene.
 *
 * This function processes animations from the provided `aiScene` and builds a mapping between armature names
 * and their corresponding animation indices. It ensures that only animations with correctly formatted names
 * are included in the mapping. The naming convention for armatures and animations is as follows:
 *
 * - Armatures must be named in the format: `X_..._armature`
 *   - `X` represents a unique identifier or descriptive name for the armature (e.g., `bottom_cylinder`, `robot`).
 *   - The `_armature` suffix is mandatory to identify the object as an armature.
 *
 * - Animations (actions) must be named in the format: `Y_..._anim`
 *   - `Y` represents a unique identifier or descriptive name for the animation action (e.g., `bottom_cylinder`,
 * `robot_idle`).
 *   - The `_anim` suffix is mandatory to identify the object as an animation action.
 *
 * - The mapping includes only those armature-animation pairs where the prefixes match:
 *   - Example: `bottom_cylinder_armature` and `bottom_cylinder_anim`
 *   - If the armature prefix does not match the animation prefix, it will not be included in the mapping.
 *
 * @param scene Pointer to the `aiScene` object containing the animations.
 * @return A mapping (`std::unordered_map<std::string, std::vector<unsigned int>>`)
 *         where the keys are armature names, and the values are vectors of animation indices.
 *         Returns an empty map if the scene is null or if no matching pairs are found.
 */
/**
 * @brief Builds a map between unique armature names and their corresponding animation indices.
 *
 * @param scene Pointer to the aiScene containing animations.
 * @return std::unordered_map<std::string, unsigned int> A map from armature names to animation indices.
 * @throws std::runtime_error If duplicate armature names are detected.
 */
std::unordered_map<std::string, unsigned int> build_armature_to_animation_map(const aiScene *scene) {
    if (!scene) {
        std::cerr << "Invalid scene pointer provided!" << std::endl;
        return {};
    }

    std::cout << "Building armature to animation map..." << std::endl;

    std::unordered_map<std::string, unsigned int> armature_to_animation_map;
    std::regex animation_name_regex(R"((.+)_armature\|(.+)_anim)");

    for (unsigned int animation_index = 0; animation_index < scene->mNumAnimations; ++animation_index) {
        aiAnimation *animation = scene->mAnimations[animation_index];
        std::string animation_name = animation->mName.C_Str();
        std::smatch match;

        std::cout << "Processing animation " << animation_index << ": " << animation_name << std::endl;

        // Check if the animation name matches the expected pattern
        if (std::regex_match(animation_name, match, animation_name_regex) && match.size() == 3) {
            std::string armature_name = match[1].str() + "_armature";
            std::string action_name = match[2].str();

            std::cout << "  Parsed armature name: " << armature_name << std::endl;
            std::cout << "  Parsed action name: " << action_name << std::endl;

            // Strip "_armature" from the armature name to compare with action_name
            std::string armature_base_name = match[1].str(); // Use the raw match group directly

            std::cout << "  Stripped armature name: " << armature_base_name << std::endl;

            // Ensure action_name starts with armature_base_name
            if (action_name.find(armature_base_name) == 0) { // Check if action_name starts with armature_base_name
                std::cout << "  Armature and action names match." << std::endl;

                // Check for duplicate armature names
                if (armature_to_animation_map.find(armature_name) != armature_to_animation_map.end()) {
                    std::cerr << "Error: Duplicate armature name detected: " << armature_name << std::endl;
                    throw std::runtime_error("Duplicate armature name detected: " + armature_name);
                }

                // Map the unique armature name to the animation index
                armature_to_animation_map[armature_name] = animation_index;
                std::cout << "  Added to map: " << armature_name << " -> " << animation_index << std::endl;
            } else {
                std::cout << "  Skipped: Armature and action names do not match." << std::endl;
            }
        } else {
            std::cout << "  Skipped: Name does not match pattern." << std::endl;
        }
    }

    std::cout << "Completed building armature to animation map. Total entries: " << armature_to_animation_map.size()
              << std::endl;

    return armature_to_animation_map;
}

/**
 * @brief Checks if the given aiNode's name matches the format X_armature.
 *
 * The function validates that the node's name ends with `_armature` and
 * has a non-empty prefix `X` before the suffix.
 *
 * @param node The aiNode to check.
 * @return True if the node's name matches the format X_armature, false otherwise.
 */
bool is_armature_node(const aiNode *node) {
    if (!node) {
        return false; // Invalid node
    }

    // Convert aiString to std::string
    std::string node_name = node->mName.C_Str();

    // Regex to match the pattern X_armature, where X is any non-empty string
    std::regex armature_regex(R"((.+)_armature)");
    return std::regex_match(node_name, armature_regex);
}

void print_all_animations(const aiScene *scene) {
    if (!scene) {
        std::cerr << "Invalid scene!" << std::endl;
        return;
    }

    // Iterate over all animations and print their names
    for (unsigned int i = 0; i < scene->mNumAnimations; ++i) {
        aiAnimation *animation = scene->mAnimations[i];
        std::cout << "Animation " << i << ": " << animation->mName.C_Str() << std::endl;
        print_ai_animation(animation);
    }
}

unsigned int find_animation_index_by_name(const aiScene *scene, const std::string &animationName) {
    // Loop through the animations to find the index of the animation with the given name
    for (unsigned int i = 0; i < scene->mNumAnimations; ++i) {
        if (scene->mAnimations[i]->mName.C_Str() == animationName) {
            return i; // Return the index if found
        }
    }
    return -1; // Return an invalid index if animation is not found
}

const aiNodeAnim *find_node_anim(const aiAnimation *pAnimation, const std::string &NodeName) {
    for (unsigned int i = 0; i < pAnimation->mNumChannels; i++) {
        const aiNodeAnim *pNodeAnim = pAnimation->mChannels[i];

        if (std::string(pNodeAnim->mNodeName.data) == NodeName) {
            return pNodeAnim;
        }
    }
    return NULL;
}

unsigned int find_idx_of_scaling_key_for_given_time(float animation_time_ticks, const aiNodeAnim *node_anim) {
    assert(node_anim->mNumScalingKeys > 0);

    for (unsigned int i = 0; i < node_anim->mNumScalingKeys - 1; i++) {
        float t = (float)node_anim->mScalingKeys[i + 1].mTime;
        if (animation_time_ticks < t) {
            return i;
        }
    }

    return 0;
}

unsigned int find_idx_of_rotation_key_for_given_time(float animation_time_ticks, const aiNodeAnim *node_anim) {
    assert(node_anim->mNumRotationKeys > 0);

    for (unsigned int i = 0; i < node_anim->mNumRotationKeys - 1; i++) {
        float t = (float)node_anim->mRotationKeys[i + 1].mTime;
        if (animation_time_ticks < t) {
            return i;
        }
    }

    return 0;
}

unsigned int find_idx_of_translation_key_for_given_time(float animation_time_ticks, const aiNodeAnim *node_anim) {
    assert(node_anim->mNumPositionKeys > 0);

    for (unsigned int i = 0; i < node_anim->mNumPositionKeys - 1; i++) {
        float t = (float)node_anim->mPositionKeys[i + 1].mTime;
        if (animation_time_ticks < t) {
            return i;
        }
    }

    return 0;
}

void RecIvpntRiggedCollector::update_animation_matrices(float animation_time_ticks) {
    rec_update_animation_matrices(animation_time_ticks, glm::mat4(1.0f), this->scene->mRootNode, this->scene, 0);
}

// Note that this data is state and contains information about the vertices of the mesh, that only need to
// be computed exactly one time, this data should get buffered into opengl one time.
std::vector<IVPNTRigged> RecIvpntRiggedCollector::parse_model_into_ivpntrs(const std::string &model_path) {
    recursion_level_counter = 0;
    const aiScene *scene = this->importer.ReadFile(model_path, aiProcess_Triangulate | aiProcess_CalcTangentSpace);
    this->scene = scene;

    if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
        std::cerr << "Error: Assimp - " << importer.GetErrorString() << std::endl;
    }

    this->directory_to_asset_being_loaded = model_path.substr(0, model_path.find_last_of("/") + 1);

    glm::mat4 root_node_transform = ai_matrix4x4_to_glm_mat4(scene->mRootNode->mTransformation);
    print_matrix(root_node_transform, "root_node_transform");
    inverse_root_node_transform = glm::inverse(root_node_transform);
    print_matrix(inverse_root_node_transform, "inverse_root_node_transform");

    print_all_animations(scene);
    armature_node_name_to_animation_index = build_armature_to_animation_map(scene);
    print_armature_to_animation_map(armature_node_name_to_animation_index);

    this->rec_process_nodes(scene->mRootNode, scene);
    return this->ivpntrs;
}

IVPNTRigged RecIvpntRiggedCollector::process_mesh_ivpntrs(aiMesh *mesh, const aiScene *scene) {
    /*std::vector<glm::vec3> vertices = process_mesh_vertex_positions(mesh, this->swap_y_and_z);*/
    std::vector<unsigned int> indices = model_loading::process_mesh_indices(mesh);
    std::vector<glm::vec3> vertices = model_loading::process_mesh_vertex_positions(mesh);
    std::vector<glm::vec3> normals = model_loading::process_mesh_normals(mesh);
    std::vector<glm::vec2> texture_coordinates = model_loading::process_mesh_texture_coordinates(mesh);
    std::vector<model_loading::TextureInfo> texture_data =
        model_loading::process_mesh_materials(mesh, scene, this->directory_to_asset_being_loaded);
    std::string main_texture = texture_data[0].path;
    std::filesystem::path fs_path = main_texture;
    // Convert to the preferred format for the operating system
    std::string texture_native_path = fs_path.make_preferred().string();

    std::vector<VertexBoneData> bone_data = this->process_mesh_vertices_bone_data(mesh);

    return {indices, vertices, normals, texture_coordinates, texture_native_path, bone_data};
};

int RecIvpntRiggedCollector::get_next_bone_id(const aiBone *pBone) {

    std::string bone_name(pBone->mName.C_Str());
    int bone_id = 0;

    if (bone_name_to_unique_index.find(bone_name) == bone_name_to_unique_index.end()) {
        bone_id = static_cast<int>(bone_name_to_unique_index.size());
        bone_name_to_unique_index[bone_name] = bone_id;
    } else {
        bone_id = bone_name_to_unique_index[bone_name];
    }

    return bone_id;
}

std::vector<VertexBoneData> RecIvpntRiggedCollector::process_mesh_vertices_bone_data(aiMesh *mesh) {
    // initialize the vector with one vertexbonedata object per vertex
    std::vector<VertexBoneData> bone_data_for_mesh(mesh->mNumVertices);
    std::cout << "working on bones of a mesh now it has: " << mesh->mNumBones << "bones" << std::endl;
    for (unsigned int i = 0; i < mesh->mNumBones; i++) {
        auto bone = mesh->mBones[i];
        std::cout << "Bone '" << bone->mName.C_Str() << "' affects " << bone->mNumWeights << " vertices" << std::endl;

        int bone_id = get_next_bone_id(bone);

        // whenever you get a new bone_id it is either reused or it is the next one
        // whenever it is the next one and not being reused, then we are looking at a new bone
        // therefore we should add to the bone_info thing
        if (bone_id == bone_unique_idx_to_info.size()) {
            BoneInfo bi(ai_matrix4x4_to_glm_mat4(bone->mOffsetMatrix));
            print_matrix(ai_matrix4x4_to_glm_mat4(bone->mOffsetMatrix), "bone offset matrix");
            bone_unique_idx_to_info.push_back(bi);
        }

        // for each bone it has a list of weights for each vertex that it affects.
        // note that this is in th eopposite order as the opengl pipeline where we are given a vertex
        // thus we require all this additionall infrastructure
        for (unsigned int j = 0; j < bone->mNumWeights; j++) { // Changed inner loop index to 'j'
            const aiVertexWeight &vw = bone->mWeights[j];
            unsigned int index_of_vertex_influenced_by_this_bone = vw.mVertexId;

            /*std::cout << ">>>> the above bone influcences a vertex with id: '"*/
            /*          << index_of_vertex_influenced_by_this_bone << "' with weight " << vw.mWeight << std::endl;*/

            // Ensure the index is within bounds
            if (index_of_vertex_influenced_by_this_bone < bone_data_for_mesh.size()) {
                bone_data_for_mesh[index_of_vertex_influenced_by_this_bone].add_bone_data(bone_id, vw.mWeight);
            } else {
                std::cerr << "Warning: Vertex index out of bounds: " << index_of_vertex_influenced_by_this_bone
                          << std::endl;
            }
        }
    }

    return bone_data_for_mesh;
}

void VertexBoneData::add_bone_data(unsigned int BoneID, float Weight) {
    for (unsigned int i = 0; i < 4; i++) {
        if (weight_value_of_this_vertex_wrt_bone[i] == 0.0) {
            indices_of_bones_that_affect_this_vertex[i] = BoneID;
            weight_value_of_this_vertex_wrt_bone[i] = Weight;
            /*std::cout << "Bone ID " << BoneID << " weight " << Weight << " stored at local index " << i <<
             * std::endl;*/
            return;
        }
    }

    bool logging = false;
    if (logging) {
        std::cout << "was about to add bone data, but we've already associated 4 weights, not adding" << std::endl;
    }
    /*assert(false); // Should never get here if we have enough space for bones, otherwise we need to increment the
     * num*/
    /*// bones count*/
}
